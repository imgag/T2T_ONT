#!/usr/bin/env python3
"""
Extract software versions for T2T-ONT pipeline tools.
Run from the project root directory:
    python workflow/scripts/get_tool_versions.py

For each conda package it tries (in order):
  1. Version pinned in the env yml file.
  2. `conda list -p <env_path> <pkg>` (fast, no subprocess activation).
  3. `conda run -p <env_path> <exe> --version` (for packages not in conda registry).
Binary tools are queried directly from the paths in config.yml.
"""

import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML required: pip install pyyaml")

PROJECT_ROOT = Path(".")
ENV_DIR      = PROJECT_ROOT / "workflow" / "env"
CONFIG_FILE  = PROJECT_ROOT / "workflow" / "config.yml"

# Additional snakemake conda cache directories to scan
SNAKEMAKE_CONDA_DIRS = [
    PROJECT_ROOT / ".snakemake" / "conda",
    PROJECT_ROOT / "analysis_other" / "cenmap" / ".snakemake" / "conda",
    Path("/mnt/storage5/users/ahgrosc1/.snakemake-conda"),
]

# ── Which env files each workflow uses ────────────────────────────────────────
WORKFLOW_ENVS = {
    "assembly":          ["verkko", "merqury", "yak", "minimap2"],
    "basecalling":       ["pod5"],
    "data_preparation":  ["filtlong", "minimap2", "py_samtools"],
    "error_correction":  [],   # uses HERRO singularity — handled in BINARY_TOOLS
    "mapping":           ["minimap2", "samtools"],
    "qc_input":          ["fastcat", "mosdepth", "R"],
    "qc_assembly":       ["merqury", "minimap2", "samtools", "bandage",
                          "gqc", "repeatmasker", "liftoff", "cenmap", "nucflag",
                          "R", "py_report"],
    "qc_extended":       ["bcftools", "mosdepth", "samtools", "whatshap", "py_report"],
    "call_variants":     ["bcftools", "whatshap", "mosdepth"],
    "methylation":       [],   # uses modkit/longphase binaries — handled in BINARY_TOOLS
    "parent_of_origin":  ["bcftools"],  # patmat binary handled separately
    "polishing":         ["medaka"],
    "dip3d":             ["bcftools", "mosdepth", "samtools", "whatshap", "py_report"],
    "porec":             ["cooler", "cooltools", "hicexplorer", "pairtools", "py_report"],
}

# Packages that are infrastructure/language runtimes, not pipeline tools
SKIP_PKGS = {
    "pip", "tabix", "pigz", "pairix",
    "pod5", "polars", "polars-lts-cpu",  # file format library, not pipeline tool
    "r-argparse", "r-optparse", "r-proto", "r-scales",
    "r-data.table", "r-viridis", "r-ggplot2", "r-tidyverse",
    "r-stringr", "bioconductor-karyoploter",
    "pandas", "numpy", "scipy", "pysam", "biopython",
    "edlib", "matplotlib", "panel", "hvplot", "typer",
    "ucsc-bedgraphtobigwig", "python",
}

# ── For unpinned packages: (executable, [args]) to run inside the env ─────────
CONDA_VERSION_CMDS = {
    "bcftools":    ("bcftools",  ["--version"]),
    "samtools":    ("samtools",  ["--version"]),
    "bandage_ng":  ("Bandage",   ["--version"]),
    "cooler":      ("cooler",    ["--version"]),
    "cooltools":   ("cooltools", ["--version"]),
    "hicexplorer": ("hicBuildMatrix", ["--version"]),
    "bedtools":    ("bedtools",  ["--version"]),
    "winnowmap":   ("winnowmap", ["--version"]),
    "liftoff":     ("liftoff",   ["--version"]),
    "nucflag":     ("nucflag",   ["--version"]),
    "r-base":      ("R",         ["--version"]),
    "fastcat":     ("fastcat",   ["--version"]),
    "filtlong":    ("filtlong",  ["--version"]),
    "medaka":      ("medaka",    ["--version"]),
    "mosdepth":    ("mosdepth",  ["--version"]),
    "whatshap":    ("whatshap",  ["--version"]),
    "pairtools":   ("pairtools", ["--version"]),
}

# ── Binary tools outside conda ────────────────────────────────────────────────
# (config_key, version_args, display_label, [workflows])
BINARY_TOOLS = [
    # (config_key,    version_args,    label,       [workflows])
    ("dorado",        ["--version"],   "dorado",         ["basecalling"]),
    (None,            None,            "herro",          ["error_correction"]),  # singularity
    ("hifiasm",       ["--version"],   "hifiasm",        ["assembly"]),
    ("mashmap",       ["--version"],   "mashmap",        ["qc_assembly"]),
    ("longdust",      ["-v"],          "longdust",       ["qc_assembly"]),
    ("falign",        None,            "falign",         ["dip3d"]),
    ("seqkit_path",   ["version"],     "seqkit",         ["qc_assembly"]),
    ("tidk_path",     ["--version"],   "tidk",           ["qc_assembly"]),
    ("covcal_path",   None,            "cov_cal",        ["qc_assembly"]),
    (None,            None,            "GFAse",          ["assembly"]),
    ("modkit",        ["--version"],   "modkit",         ["methylation"]),
    ("longphase",     ["--version"],   "longphase",      ["methylation"]),
    ("patmat",        None,            "patmat",         ["parent_of_origin"]),
    (None,            ["--version"],   "alfred",         ["qc_assembly"]),
]
GFASE_PATH   = "bin/GFAse/build/phase_contacts_with_monte_carlo"
HERRO_LABEL  = "herro-v1 (singularity)"
ALFRED_PATH  = "bin/alfred"

CONTAINER_TOOLS = [
    ("hmm_flagger", ["assembly"]),
]

# Matches semver-like strings (require a dot separator to avoid matching bare dates)
VERSION_RE = re.compile(r"(\d+\.\d[\d\.]*(?:-r\d+)?)")


def run(cmd, timeout=10):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return (r.stdout + r.stderr).strip()
    except Exception:
        return None


def first_version(text):
    if not text:
        return None
    # Strip log-prefix lines (e.g. "[2026-07-24 13:59:36.882] [info] ...")
    lines = [l for l in text.splitlines() if not re.match(r"^\s*\[?\d{4}-\d{2}-\d{2}", l)]
    cleaned = "\n".join(lines) if lines else text
    m = VERSION_RE.search(cleaned)
    return m.group(1) if m else cleaned.split("\n")[0][:40]


# ── Build env name -> installed path map ──────────────────────────────────────

def _build_env_map():
    """
    Return {env_name: absolute_path} by scanning:
      1. All snakemake conda cache dirs (match via name: field in yaml).
      2. Named envs from `mamba env list` / `conda env list`.
    """
    env_map = {}

    # 1. Snakemake hash envs
    for conda_dir in SNAKEMAKE_CONDA_DIRS:
        if not conda_dir.exists():
            continue
        for yaml_file in conda_dir.glob("*.yaml"):
            env_dir = yaml_file.parent / yaml_file.stem  # strips .yaml
            if not env_dir.exists():
                continue
            try:
                with open(yaml_file) as f:
                    data = yaml.safe_load(f) or {}
                name = data.get("name", "")
                if name and name not in env_map:
                    env_map[name] = str(env_dir)
            except Exception:
                continue

    # 2. Named envs from mamba/conda env list
    for cmd in (["mamba", "env", "list"], ["conda", "env", "list"]):
        out = run(cmd, timeout=15)
        if not out:
            continue
        for line in out.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2 and not parts[0].startswith("/"):
                name, path = parts[0], parts[-1]
                if name not in env_map:
                    env_map[name] = path
        break  # stop after first successful command

    return env_map


_ENV_MAP = None


def env_path(name):
    global _ENV_MAP
    if _ENV_MAP is None:
        _ENV_MAP = _build_env_map()
    return _ENV_MAP.get(name)


# ── Version resolution ────────────────────────────────────────────────────────

def _conda_list_at(path, pkg_name):
    """Return version string from `conda list -p <path> <pkg>`, or None."""
    if not path:
        return None
    out = run(["conda", "list", "-p", path, pkg_name, "--json"], timeout=15)
    if not out:
        return None
    try:
        for p in json.loads(out):
            if p.get("name", "").lower() == pkg_name.lower():
                return p.get("version")
    except Exception:
        pass
    return None


def conda_list_version(env_name, pkg_name):
    """Query installed version, trying the yml env name then the package name as env."""
    # Primary: env named after the yml file
    v = _conda_list_at(env_path(env_name), pkg_name)
    if v:
        return v
    # Fallback: env whose name matches the package (e.g. 'bandage_ng' env for 'bandage_ng' pkg)
    if pkg_name != env_name:
        v = _conda_list_at(env_path(pkg_name), pkg_name)
    return v


def conda_run_version(env_name, pkg_name):
    """Run the tool inside the env to get its version string."""
    if pkg_name not in CONDA_VERSION_CMDS:
        return None
    exe, args = CONDA_VERSION_CMDS[pkg_name]
    # Try yml env name, then package name as env
    for name in ([env_name] if env_name == pkg_name else [env_name, pkg_name]):
        path = env_path(name)
        if not path:
            continue
        out = run(["conda", "run", "-p", path, exe] + args, timeout=15)
        v = first_version(out) if out else None
        if v:
            return v
    return None


def resolve_version(env_name, pkg_name, pinned):
    if pinned:
        return pinned
    v = conda_list_version(env_name, pkg_name)
    if v:
        return v
    v = conda_run_version(env_name, pkg_name)
    return v or "UNKNOWN"


# ── Parse env yml ─────────────────────────────────────────────────────────────

def parse_env_yml(env_name):
    path = ENV_DIR / f"{env_name}.yml"
    if not path.exists():
        return {}
    with open(path) as f:
        data = yaml.safe_load(f)
    result = {}
    for dep in data.get("dependencies", []):
        if isinstance(dep, str):
            m = re.match(r"^([A-Za-z0-9_\-\.]+)\s*[=><]+\s*([^\s]+)", dep)
            if m:
                result[m.group(1).lower()] = m.group(2).lstrip("=")
            else:
                result[dep.strip().lower()] = None
        elif isinstance(dep, dict):
            for pip_dep in dep.get("pip", []):
                m = re.match(r"^([A-Za-z0-9_\-\.]+)[=><]+\s*([^\s]+)", pip_dep)
                if m:
                    result[m.group(1).lower()] = m.group(2).lstrip("=")
                else:
                    result[pip_dep.strip().lower()] = None
    return result


# ── Collectors ────────────────────────────────────────────────────────────────

def collect_conda_tools():
    records = []
    cache = {}
    for workflow, env_names in WORKFLOW_ENVS.items():
        for env_name in env_names:
            for pkg, pinned in parse_env_yml(env_name).items():
                if pkg in SKIP_PKGS:
                    continue
                key = (env_name, pkg)
                if key not in cache:
                    cache[key] = resolve_version(env_name, pkg, pinned)
                records.append({"workflow": workflow, "tool": pkg,
                                "version": cache[key], "source": f"conda:{env_name}"})
    return records


def collect_binary_tools(config):
    records = []
    for cfg_key, ver_args, label, workflows in BINARY_TOOLS:
        if label == "GFAse":
            version, path = "dev (no release tag)", GFASE_PATH
        elif label == "herro":
            # Version encoded in directory name; singularity image in bin/herro.sif
            version = config.get("herro_model", "bin/herro-v1").split("/")[-1]
            path = "bin/herro.sif (singularity)"
        elif label == "alfred":
            path = ALFRED_PATH
            out = run([str(PROJECT_ROOT / path), "--version"])
            version = first_version(out) if out else "UNKNOWN"
        elif label == "patmat":
            path = config.get("patmat", "bin/PatMat/patmat")
            # version from pyproject.toml next to the binary
            pyproject = PROJECT_ROOT / Path(path).parent / "pyproject.toml"
            v = None
            if pyproject.exists():
                import re as _re
                m = _re.search(r'^version\s*=\s*"([^"]+)"', pyproject.read_text(), _re.M)
                v = m.group(1) if m else None
            version = v or "UNKNOWN"
        elif cfg_key is None:
            continue
        else:
            path = config.get(cfg_key, cfg_key)
            if ver_args is None:
                version = "UNKNOWN (no --version)"
            else:
                out = run([str(PROJECT_ROOT / path)] + ver_args)
                version = first_version(out) if out else "UNKNOWN"
        for wf in workflows:
            records.append({"workflow": wf, "tool": label,
                            "version": version, "source": f"binary:{path}"})
    return records


def collect_container_tools(config):
    records = []
    for cfg_key, workflows in CONTAINER_TOOLS:
        val = config.get(cfg_key, "")
        m = re.search(r":v?([\d\.]+)", val)
        version = m.group(1) if m else val
        for wf in workflows:
            records.append({"workflow": wf, "tool": cfg_key,
                            "version": version, "source": f"container:{val}"})
    return records


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)

    records = (collect_conda_tools()
               + collect_binary_tools(config)
               + collect_container_tools(config))

    # Deduplicate: merge workflows for identical (tool, source) pairs
    tool_info = defaultdict(lambda: {"workflows": set(), "version": "", "source": ""})
    for r in records:
        key = (r["tool"], r["source"])
        tool_info[key]["version"]  = r["version"]
        tool_info[key]["source"]   = r["source"]
        tool_info[key]["workflows"].add(r["workflow"])

    print("\t".join(["tool", "version", "workflows", "source"]))
    for (tool, source), info in sorted(tool_info.items(), key=lambda x: x[0][0].lower()):
        wf = ",".join(sorted(info["workflows"]))
        print(f"{tool}\t{info['version']}\t{wf}\t{source}")


if __name__ == "__main__":
    main()
