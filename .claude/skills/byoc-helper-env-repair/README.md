# byoc helper-env repair

Rebuilds the **local** conda env that a Claude Science compute provider's
`compute_provider` kernel runs inside, for the case where the platform's own
provisioning of that env failed.

Nothing here touches your cloud account. The env holds the provider SDK
(`modal`, `httpx`, …) on your machine; when it's missing, the kernel exits
before any remote call happens.

## The failure it fixes

A `compute_provider` cell dies at startup with:

```
ByocError: ('network_bridge_down',
  'aiohttp trust_env patch failed: ModuleNotFoundError("No module named 'modal'")
   — blob I/O would bypass the proxy')
```

`network_bridge_down` is misleading: it's a **missing package**, not a network
problem. The provider patches `modal._utils.http_utils._http_client_with_tls`
so blob I/O honors the proxy, and the patch can't import its target. The daemon
logs the same fact as `the byoc helper's python was not found
(~/.claude-science/conda/envs/compute-provider-modal/bin/python)`.

## Why provisioning failed (on this machine)

The bundled micromamba wedges on `create` before printing any solver output and
holds the platform's global **env-mutation lock**, so every queued env
mutation — including this one — starves:

```
[conda] MCP env 'claude-science-mcp': waited 620.2s for the env mutation lock
        (holder phase when the wait began: micromamba create)
```

A manual `micromamba create --dry-run` hangs identically, which is what rules
out the app as the culprit. `pgrep -fl micromamba` names the stuck holder.

## What the script does

`provision_helper_env.sh [provider-id] [runtime-dir]` — default provider
`modal`; `infer` also works.

1. Reads the authoritative spec from the installed runtime's
   `skills/*/provider.json` → `helperEnv` (env name, python version, pip
   payload), so it stays correct across app versions instead of hardcoding.
2. Reports leftover `micromamba` processes and clears any partial env.
3. **Path A**: `micromamba create --offline <pyspec> pip` — resolves only from
   `conda/pkgs`, which already holds the platform's python and pip, so the
   repodata fetch it stalls on never happens. Wrapped in a 120 s watchdog
   (macOS has no `timeout`).
4. **Path B** (if A hangs or errors): `python -m venv` off the platform's
   existing conda Python 3.11. The daemon only requires `<prefix>/bin/python`
   to exist, so a venv satisfies it.
5. Installs the payload from the provider's hash-pinned `requirements.lock`
   (`pip install --require-hashes`), giving byte-identical packages to a
   successful platform provision; falls back to the `provider.json` pip specs
   if no lock ships.
6. Stamps `.operon_metadata.json` like the other platform-managed envs, then
   verifies by importing the payload's first module — and, for Modal, asserting
   the proxy-patch target is callable.

## Running it

```bash
bash provision_helper_env.sh 2>&1 | tee ~/helper_env_provision.log
```

Run it in **Terminal**, not through the agent: the env name is reserved
(`manage_environments` refuses it) and the agent sandbox cannot write under
`~/.claude-science/conda`.

**Do not restart the app afterward.** Each `compute_provider` cell spawns a
fresh kernel, so the new env is picked up on the next cell.

## Verifying

Expected tail:

```
== interpreter: Python 3.11.16 at ~/.claude-science/conda/envs/compute-provider-modal/bin/python
   modal 1.5.1
   proxy patch target present: modal._utils.http_utils._http_client_with_tls
```

Then, agent-side: `compute_provider_config()` returns the bound Modal
Environment / app / egress mode, `app.app_id` confirms the credential against
Modal's control plane, and `list_envs()` lists the bundled GPU env definitions.

If path A times out at 120 s, that's worth reporting as a platform bug — the
hang is inside micromamba, independent of the app.

## Also available as a skill

`skill("byoc-helper-env-repair")` bundles this script plus two read-only
helpers:

- `byoc_helper_env_status(provider, runtime_dir=None)` — env name, prefix,
  whether the interpreter exists, python version, and whether the boot import
  succeeds. Without `runtime_dir` it falls back to the
  `compute-provider-<provider>` naming convention, because the agent sandbox
  can't list `~/.claude-science/runtime`.
- `byoc_helper_env_repair_command(provider)` — the exact Terminal command to
  hand to the user.
