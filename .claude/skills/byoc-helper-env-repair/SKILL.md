---
name: byoc-helper-env-repair
description: "Repair a byoc compute provider's local helper conda env when the compute_provider kernel will not start — ByocError network_bridge_down, ModuleNotFoundError No module named modal, or the ops.compute reaper reporting the byoc helper's python was not found. Covers reading the provider.json helperEnv spec, rebuilding the env offline when the bundled micromamba hangs on create and holds the env-mutation lock, installing the hash-pinned requirements.lock payload, and verifying the kernel boots. Use when compute_provider cells die at startup, when Modal or infer provisioning failed, or when conda env installs stall behind a wedged micromamba."
---

# Repairing a byoc provider's local helper env

Every byoc compute provider runs its `compute_provider` kernel inside a
**local** conda env that the platform provisions from the provider's own
spec. Nothing about the user's cloud account is involved. When that env is
missing or incomplete, the kernel dies at boot and no amount of
credential-checking or re-approving helps.

## Recognizing it

The kernel exits 1 before your code runs. For Modal the traceback ends:

```
operon_compute_provider.ByocError: ('network_bridge_down',
  'aiohttp trust_env patch failed: ModuleNotFoundError("No module named \'modal\'")
   — blob I/O would bypass the proxy')
```

`network_bridge_down` here is a **missing package**, not a network problem:
the provider patches `modal._utils.http_utils._http_client_with_tls` so blob
I/O honors the proxy, and that patch cannot import its target. The daemon
says the same thing in `~/.claude-science/logs/server-*.log`:

```
[ops.compute] periodic reaper modal pass failed: the byoc helper's python was not found
 (~/.claude-science/conda/envs/compute-provider-modal/bin/python) — ... not yet provisioned
```

## Confirm the spec, then the disk

The provider declares what the env must contain in
`~/.claude-science/runtime/<ver>/skills/<provider-skill>/provider.json` under
`helperEnv` — env name, python version, pip payload. Two examples shipped in
0.1.51: `remote-compute-modal` → `compute-provider-modal`, `python=3.11` +
`pip`, then `modal==1.5.1`, `python-socks[asyncio]==2.8.1`,
`aiohttp-socks==0.11.0`; `using-model-endpoint` → `compute-provider-infer`,
`httpx==0.28.1`. The sibling `requirements.lock` is the hash-pinned closure
the platform itself installs.

`byoc_helper_env_status("modal")` (from this skill's `kernel.py`) reads that
spec, checks whether `<conda root>/envs/<name>/bin/python` exists, and tries
the import the kernel does at boot. Read-only, safe to call any time.

The agent sandbox cannot **list** `~/.claude-science/runtime` (it can read
inside a version dir once named), so the `provider.json` glob comes back empty
there and the helper falls back to the `compute-provider-<provider>` naming
convention — enough for the exists/import check. Pass
`runtime_dir="~/.claude-science/runtime/<ver>"` (the version is printed in the
kernel traceback) to get the pip payload and lock path too. The script takes
the same value as its optional second argument.

## Why provisioning fails, and the fix

The common cause is not the provider: the bundled micromamba wedges on
`create` before emitting solver output and **holds the platform's global
env-mutation lock**, so every queued mutation — including this env — starves.
The signature in `server-*.log`:

```
[conda] MCP env 'claude-science-mcp': waited 620.2s for the env mutation lock
        (previous holder's last phase: entering mutation;
         holder phase when the wait began: micromamba create)
```

A manual `micromamba create --dry-run` hangs the same way, which confirms the
hang is in micromamba rather than the app. `pgrep -fl micromamba` names the
holder; the user quitting the app (or killing that PID) releases it and also
unblocks whatever other env work was queued.

The repair avoids the network path micromamba stalls on:

1. `micromamba create -r <root> -p <prefix> --offline <pyspec> pip` — resolves
   from `conda/pkgs`, which already holds the platform's own python and pip,
   so no repodata fetch.
2. Fallback if that hangs or errors: `<root>/envs/python/bin/python3.11 -m venv
   <prefix>`. The daemon only requires `<prefix>/bin/python` to exist, so a
   venv satisfies it just as well.
3. `pip install --require-hashes -r <provider skill>/requirements.lock`, giving
   byte-identical packages to a successful platform provision.
4. Stamp `.operon_metadata.json` (`{language, op_log:[...]}`) like the other
   platform-managed envs.
5. Verify by importing the payload's first module and, for Modal, asserting
   `modal._utils.http_utils._http_client_with_tls` is callable.

`scripts/provision_helper_env.sh` does all five, parameterized by provider id
(default `modal`), reading the spec from the installed runtime so it stays
correct across versions. It wraps step 1 in a 120 s watchdog because macOS
ships no `timeout`.

## Two constraints that shape the workflow

**You cannot run the fix yourself.** The helper env name is reserved
(`manage_environments` refuses it: "reserved environment name
(system-owned)") and the agent sandbox cannot write under
`~/.claude-science/conda` (`Operation not permitted`), so the script runs in
the user's Terminal. Hand them the command from
`byoc_helper_env_repair_command(provider)` and ask them to paste the output.

**No app restart.** Each `compute_provider` cell spawns a fresh kernel, so the
new env is picked up on the next cell. Say this explicitly — a restart wastes
a cycle, and if the env vanishes across one, suspect a startup sweep and
rebuild without restarting.

## After it comes up

Confirm reach before building anything: `compute_provider_config()` for the
bound Environment / app / `egress_mode`, `app.app_id` for a real control-plane
round trip on the user's credential, and `list_envs()` for the bundled image
definitions. Then append a short gotcha note to `compute_details` for that
provider so the next session doesn't rediscover the hang.
