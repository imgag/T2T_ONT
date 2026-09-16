# T2T_ONT (diagnostic branch)

Snakemake pipeline for telomere-to-telomere (T2T) genome assembly and QC from
ONT long-read sequencing, scoped to **diagnostic samples only** (`T2TP*`
sample IDs). This is the `diagnostic` branch of the `T2T_ONT` GitHub repo
(`imgag/T2T_ONT`).

## Relationship to the sister project

This is a sister project to `../ahthapp1_T2T_ONT` (same GitHub repo,
`imgag/T2T_ONT`, but tracking `origin` = `imgag/T2T_ONT_dev` there). That
project is the main research/development repo (assembly QC, ancestry,
Pore-C, 3D structure prediction, publication figures, etc.) and sees far more
active development. This repo intentionally only carries the subset of the
pipeline needed to process diagnostic samples, plus their pipeline-level
(non-sample-specific) fixes as they land upstream — it does not track every
research feature added there.

When porting changes from `../ahthapp1_T2T_ONT`:
- Prefer narrow, pipeline-level changes (basecalling, dorado, error
  correction, core QC) over research-feature additions (ancestry, 3D
  structure, Pore-C-specific analyses) unless a diagnostic sample actually
  needs them.
- Never overwrite `data/datasets.yml`, `data/finished_samples.yml`,
  `assembly/assemblies.yml`, or anything under `data/raw`/`data/basecalled` —
  these are this repo's own sample tracking and must be preserved across any
  merge.
- Some config values intentionally diverge from the sister repo (e.g. a
  `dorado_model.mod` version that isn't downloaded here) — check that a
  referenced model/binary path actually exists under `bin/` or
  `data/dorado_models/` before copying a config value over.

## Structure

- `workflow/Snakefile` — entry point; loads `workflow/config.yml`,
  `data/datasets.yml`, `assembly/assemblies.yml`, `data/finished_samples.yml`,
  then includes all rule files from `workflow/rules/`.
- `workflow/rules/` — one `.smk` file per pipeline stage: `basecalling.smk`
  (dorado basecaller/duplex), `error_correction.smk` (dorado trim/correct +
  herro), `data_preparation.smk` (dataset/input-path resolution helpers),
  `mapping.smk`, `assembly.smk` (verkko/hifiasm), `polishing.smk`,
  `qc_input.smk`, `qc_assembly.smk`, `qc_extended.smk`, `call_variants.smk`,
  `dip3d.smk`, `porec.smk`, `ancestry.smk`, `_functions.smk` (shared helper
  functions), `handle_gpu.smk` (GPU slot allocation for dorado jobs).
- `workflow/env/` — one conda env per tool.
- `workflow/scripts/` — R/Python scripts invoked from rules (QC plots, dorado
  summary stats).
- `data/datasets.yml` — maps each diagnostic sample (`T2TP*`) to its raw
  run-folder inputs per data type (`UL`, `HQ_herro`, `POREC`, ...). This is
  sample tracking, not pipeline code — keep it as-is unless explicitly asked
  to add/update a sample.
- `data/finished_samples.yml`, `assembly/assemblies.yml` — further sample
  state tracking, same rule applies.
- `bin/` — local tool binaries (dorado, hifiasm, etc.), gitignored. Some
  config entries point here with a local version, others point across to
  `../ahthapp1_T2T_ONT/bin/...` for tools not duplicated locally — follow
  whichever convention the existing config key already uses.
- `doc/` — analysis documentation and QC writeups.

## Working with this repo

- GPU-bound dorado rules (`dorado`, `dorado_duplex`, `dorado_correct_inference`)
  acquire a GPU slot via `handle_gpu.smk`'s `get_gpu_id()` and use
  `resources: queue=config['gpu_queues']` — keep new GPU rules consistent
  with this pattern rather than hardcoding a queue name.
- Recent dorado versions require
  `LD_LIBRARY_PATH="$(dirname <dorado_bin>)/../lib:$LD_LIBRARY_PATH"` before
  invocation — this is already wired into the basecalling/error-correction
  shell blocks; keep it when touching those rules.
- `bin/` and most raw/intermediate data directories are gitignored (see
  `.gitignore`) — large binaries and sample data are not meant to be
  committed.
