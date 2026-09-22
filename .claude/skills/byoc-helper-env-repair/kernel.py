import glob
import json
import os


def byoc_helper_env_status(provider="modal", runtime_dir=None):
    """Read-only status of a byoc provider's LOCAL helper conda env.

    Returns env_name / prefix / exists / python_version / import_ok, plus the
    provider.json and requirements.lock paths when the spec is readable.

    Path note: the agent sandbox cannot LIST ~/.claude-science/runtime (though
    it can read inside a known version dir), so the provider.json glob fails
    there. Pass runtime_dir="~/.claude-science/runtime/<ver>" when you know the
    version (kernel tracebacks print it); otherwise the env name falls back to
    the convention "compute-provider-<provider>".
    """
    import subprocess
    cs = os.path.expanduser("~/.claude-science")
    if runtime_dir:
        pattern = os.path.join(os.path.expanduser(runtime_dir), "skills", "*", "provider.json")
    else:
        pattern = os.path.join(cs, "runtime", "*", "skills", "*", "provider.json")
    hits = []
    for path in sorted(glob.glob(pattern)):
        try:
            spec = json.load(open(path))
        except Exception:
            continue
        if spec.get("id") == provider:
            hits.append((path, spec))
    out = {"provider": provider}
    if hits:
        path, spec = hits[-1]
        helper = spec.get("helperEnv") or {}
        out["provider_json"] = path
        out["lock"] = os.path.join(os.path.dirname(path), "requirements.lock")
        out["pip"] = helper.get("pip", [])
        out["env_name"] = helper.get("name") or ("compute-provider-" + provider)
        out["spec_source"] = "provider.json"
    else:
        out["env_name"] = "compute-provider-" + provider
        out["pip"] = []
        out["spec_source"] = "convention (provider.json not readable from here)"
    prefix = os.path.join(cs, "conda", "envs", out["env_name"])
    py = os.path.join(prefix, "bin", "python")
    out["prefix"] = prefix
    out["exists"] = os.path.exists(py)
    if not out["exists"]:
        return out
    if out["pip"]:
        first = out["pip"][0]
        for sep in ("==", ">=", "<=", "~=", "["):
            first = first.split(sep)[0]
        module = first.strip().replace("-", "_")
    else:
        module = {"modal": "modal", "infer": "httpx"}.get(provider, provider)
    out["import_module"] = module
    out["python_version"] = subprocess.run(
        [py, "-V"], capture_output=True, text=True).stdout.strip()
    probe = subprocess.run(
        [py, "-c", "import importlib,sys;m=importlib.import_module(sys.argv[1]);"
                   "print(getattr(m,'__version__',''))", module],
        capture_output=True, text=True)
    out["import_ok"] = probe.returncode == 0
    out["import_detail"] = (probe.stdout or probe.stderr).strip()[-300:]
    return out


def byoc_helper_env_repair_command(provider="modal", skill_dir=None):
    """The Terminal command for the USER to run (the sandbox cannot write
    under ~/.claude-science/conda). Returns a string to hand over verbatim."""
    import sys
    if skill_dir is None:
        skill_dir = os.path.dirname(sys._getframe().f_code.co_filename)
    script = os.path.join(skill_dir, "scripts", "provision_helper_env.sh")
    if not skill_dir:
        script = "provision_helper_env.sh"
    return "bash %s %s 2>&1 | tee ~/helper_env_provision.log" % (script, provider)
