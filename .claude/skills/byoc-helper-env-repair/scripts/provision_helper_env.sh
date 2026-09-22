#!/bin/bash
# Rebuild a byoc compute-provider's local helper conda env by hand, without
# letting micromamba touch the network.
#
#   bash provision_helper_env.sh [provider-id] [runtime-dir]
#     provider-id  default: modal  (also: infer)
#     runtime-dir  optional, e.g. ~/.claude-science/runtime/0.1.51-release —
#                  only needed if the runtime glob can't see the version dirs
#
# Run in Terminal, not in the agent sandbox — the sandbox cannot write under
# ~/.claude-science/conda. No app restart is needed afterward: each
# compute_provider cell spawns a fresh kernel.
#
# Reads the provider's own helperEnv declaration out of
# runtime/<ver>/skills/*/provider.json, so the env name, python version and
# pip payload always match what the installed runtime expects.

set -euo pipefail

PROVIDER="${1:-modal}"
RUNTIME_DIR="${2:-}"
CS="$HOME/.claude-science"
ROOT="$CS/conda"
MM="$ROOT/bin/micromamba"

BASEPY="$ROOT/envs/python/bin/python3.11"
[ -x "$BASEPY" ] || BASEPY="$ROOT/envs/python/bin/python3"
[ -x "$BASEPY" ] || BASEPY="$(command -v python3 || true)"
[ -n "$BASEPY" ] && [ -x "$BASEPY" ] || { echo "FAIL: no python3 to read provider.json with"; exit 1; }

# ---- read the provider's helperEnv spec -------------------------------------
TMPPY="$(mktemp -t byocprov)"
cat > "$TMPPY" <<'PY'
import glob, json, os, sys
cs, prov = sys.argv[1], sys.argv[2]
rt = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else None
pattern = (os.path.join(os.path.expanduser(rt), "skills", "*", "provider.json") if rt
           else os.path.join(cs, "runtime", "*", "skills", "*", "provider.json"))
pick = None
for p in sorted(glob.glob(pattern)):
    try:
        d = json.load(open(p))
    except Exception:
        continue
    if d.get("id") == prov:
        pick = (p, d)          # last match wins == newest runtime
if pick is None:
    print("NONE NONE NONE NONE")
    raise SystemExit(0)
path, d = pick
he = d.get("helperEnv") or {}
pyspec = next((x for x in he.get("packages", []) if x.startswith("python")), "python=3.11")
pips = ",".join(he.get("pip", []) or ["NONE"])
print(he.get("name", "NONE"), pyspec, os.path.dirname(path), pips)
PY
read -r ENV_NAME PYSPEC PROVDIR PIPSPECS <<<"$("$BASEPY" "$TMPPY" "$CS" "$PROVIDER" "$RUNTIME_DIR")"
rm -f "$TMPPY"

[ "$ENV_NAME" != "NONE" ] || { echo "FAIL: no provider.json with id='$PROVIDER' under $CS/runtime/*/skills/"; exit 1; }

PREFIX="$ROOT/envs/$ENV_NAME"
LOCK="$PROVDIR/requirements.lock"

echo "== spec (from $PROVDIR/provider.json)"
echo "   provider : $PROVIDER"
echo "   env name : $ENV_NAME"
echo "   python   : $PYSPEC"
echo "   pip      : ${PIPSPECS//,/ }"
echo "   lock     : $([ -f "$LOCK" ] && echo "$LOCK" || echo "(none — will install the pip specs directly)")"
echo "   prefix   : $PREFIX"

export MAMBA_ROOT_PREFIX="$ROOT"

# run_to <seconds> <cmd...> — macOS ships no coreutils `timeout`.
run_to() {
  local t="$1"; shift
  "$@" & local p=$!
  local i=0
  while kill -0 "$p" 2>/dev/null; do
    if [ "$i" -ge "$t" ]; then
      echo "   !! still running after ${t}s — terminating (this is the hang)"
      kill -TERM "$p" 2>/dev/null; sleep 2; kill -KILL "$p" 2>/dev/null
      wait "$p" 2>/dev/null; return 124
    fi
    sleep 1; i=$((i+1))
  done
  wait "$p"
}

echo "== leftover micromamba processes (a wedged one holds the env-mutation lock)"
pgrep -fl micromamba || echo "   none"

echo "== clearing any partial env"
if [ -e "$PREFIX" ]; then rm -rf "$PREFIX"; echo "   removed stale $PREFIX"; else echo "   nothing to clear"; fi

echo "== path A: micromamba create --offline (package cache only, no repodata fetch)"
if [ -x "$MM" ] && run_to 120 "$MM" create -y -r "$ROOT" -p "$PREFIX" --offline "$PYSPEC" pip; then
  echo "   path A ok"
else
  echo "   path A unavailable, failed, or hung — falling back"
  rm -rf "$PREFIX"
  echo "== path B: venv off the platform's existing conda python"
  echo "   base: $BASEPY ($("$BASEPY" -V 2>&1))"
  "$BASEPY" -m venv "$PREFIX"
  "$PREFIX/bin/python" -m pip install --no-input --disable-pip-version-check -q -U pip \
    || echo "   (pip self-upgrade skipped)"
fi

echo "== interpreter: $("$PREFIX/bin/python" -V 2>&1) at $PREFIX/bin/python"

echo "== installing the provider payload"
if [ -f "$LOCK" ]; then
  "$PREFIX/bin/python" -m pip install --no-input --disable-pip-version-check \
    --require-hashes -r "$LOCK"
else
  IFS=',' read -r -a SPECS <<<"$PIPSPECS"
  "$PREFIX/bin/python" -m pip install --no-input --disable-pip-version-check "${SPECS[@]}"
fi

echo "== stamping env metadata (matches the other platform-managed envs)"
"$PREFIX/bin/python" - "$PREFIX" "$PYSPEC" <<'PY'
import datetime, json, os, sys
prefix, pyspec = sys.argv[1], sys.argv[2]
meta = {"language": "python", "op_log": [{
    "timestamp": datetime.datetime.now(datetime.timezone.utc)
                 .isoformat(timespec="milliseconds").replace("+00:00", "Z"),
    "operation": "create",
    "packages": [pyspec, "pip"],
    "result": "success",
    "python_version": pyspec.split("=")[-1],
    "note": "hand-provisioned: offline micromamba create or venv; payload from provider spec",
}]}
path = os.path.join(prefix, ".operon_metadata.json")
json.dump(meta, open(path, "w"), indent=2)
print("   wrote", path)
PY

echo "== verifying the import the kernel does at boot"
IMPORT_NAME="$("$BASEPY" -c "
import sys
s = sys.argv[1].split(',')[0]
for sep in ('==', '>=', '<=', '~=', '['):
    s = s.split(sep)[0]
print(s.strip().replace('-', '_'))
" "$PIPSPECS")"
"$PREFIX/bin/python" - "$IMPORT_NAME" <<'PY'
import importlib, sys
name = sys.argv[1]
m = importlib.import_module(name)
print("  ", name, getattr(m, "__version__", "(no __version__)"))
if name == "modal":
    from modal._utils import http_utils
    assert callable(http_utils._http_client_with_tls), "proxy patch target missing"
    print("   proxy patch target present: modal._utils.http_utils._http_client_with_tls")
PY

echo
echo "== done. Do NOT restart the app — tell the agent, and it will retry the"
echo "   compute_provider kernel immediately (a fresh kernel spawns per cell)."
