# T2T-ONT Workflow

Snakemake workflow for diploid telomere-to-telomere genome assembly from Oxford Nanopore Technology (ONT) reads.

## Overview

The workflow produces haplotype-resolved, near-T2T human genome assemblies from:
- Ultra-long (UL) ONT reads
- HERRO-corrected reads (generated from simplex reads)
- Pore-C reads (for scaffolding and phasing)

It also produces phased methylation maps and haplotype-resolved 3D chromatin contact matrices from the same sequencing data.

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 8.0
- [Conda](https://docs.conda.io/) or [Mamba](https://github.com/mamba-org/mamba) (for environment management)
- GPU cluster with SLURM (for HERRO error correction and basecalling)

## Configuration

### 1. Sample configuration (`data/datasets.yml`)

Defines input data paths for each sample. Each entry specifies the data type and file paths:

```yaml
SAMPLE_NAME:
  UL:                          # Ultra-long reads (folder for basecalling or .fastq.gz)
    - "data/raw/sample/run_001"
  POREC:                       # Pore-C reads (folder or .fastq.gz)
    - "data/raw/sample/porec_001"
  HQ_herro:                    # Simplex reads for HERRO correction (.fastq.gz list)
    - "data/raw/sample/simplex_001.fastq.gz"
```

If the input path is a directory, the workflow assumes it contains raw pod5/fast5 files and runs Dorado basecalling automatically.

### 2. Assembly configuration (`assembly/assemblies.yml`)

Defines which assemblies to produce and links them to dataset entries:

```yaml
SAMPLE_NAME:
  dataset: SAMPLE_NAME         # Key in datasets.yml
  trio_phasing: false          # Enable trio-based phasing
  apk_polishing: false         # Enable APK polishing
```

### 3. Pipeline parameters (`workflow/config.yml`)

Key parameters:

| Parameter | Description |
|:----------|:------------|
| `ref` | Path to T2T-CHM13v2.0 reference FASTA |
| `dorado_model` | Dorado basecalling model paths |
| `min_length.UL` | Minimum read length for ultra-long reads (default: 80,000 bp) |
| `min_mean_q` | Minimum mean read quality (default: 9) |
| `genome_length` | Expected genome size in bp (default: 3,200,000,000) |

## Running the workflow

### Dry run (check which jobs would be executed)

```bash
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yml \
    --use-conda \
    --dryrun \
    <target>
```

### Assembly (main pipeline)

```bash
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yml \
    --use-conda \
    --cores 32 \
    <target>
```

On a SLURM cluster:

```bash
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yml \
    --use-conda \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time={resources.runtime}" \
    --jobs 50 \
    <target>
```

### Main workflow targets

| Target | Description |
|:-------|:------------|
| `assembly` | Run Verkko assembly for all configured samples |
| `qc_assembly` | Run full assembly QC (paftools, Merqury, NucFlag, Flagger) |
| `qc_input` | Run sequencing read QC |
| `polishing` | Run Medaka APK polishing |
| `call_variants` | Call small variants and SVs against CHM13 |
| `methylation` | Generate phased methylation maps |
| `porec` | Process Pore-C data and generate contact matrices |
| `ancestry` | Run population genetics analysis |

## Workflow modules

| Module | Description |
|:-------|:------------|
| `basecalling.smk` | Dorado basecalling from raw pod5/fast5 |
| `error_correction.smk` | HERRO correction of simplex reads |
| `assembly.smk` | Verkko and hifiasm assembly |
| `polishing.smk` | Medaka APK/ULK polishing |
| `qc_input.smk` | Read-level QC and statistics |
| `qc_assembly.smk` | Assembly QC (paftools, Merqury, asmgene) |
| `qc_extended.smk` | Extended QC (NucFlag, Flagger, RepeatMasker) |
| `call_variants.smk` | Variant calling (dipcall, hapdiff) |
| `methylation.smk` | Phased CpG methylation (modkit) |
| `porec.smk` | Pore-C contact matrix generation |
| `dip3d.smk` | Haplotype-resolved chromatin contacts (Dip3D) |
| `ancestry.smk` | Population genetics and ancestry inference |

## Environments

Conda environments for each step are defined in `workflow/env/`. Environments are created automatically by Snakemake on first use with `--use-conda`.
