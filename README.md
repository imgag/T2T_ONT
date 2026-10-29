# T2T-ONT: Diploid Telomere-to-Telomere Assembly from ONT Reads


This repository contains the analysis code and documentation for producing haplotype-resolved, near-T2T human genome assemblies using only Oxford Nanopore Technology (ONT) sequencing. Three types of ONT data are combined: ultra-long reads, HERRO-corrected reads, and Pore-C reads. The same sequencing data also yields phased CpG methylation maps and haplotype-resolved 3D chromatin contact matrices. It is the companion repository for this publication:

**Single-Platform Nanopore Sequencing Enables Diploid Telomere-to-Telomere Genome Assembly and Haplotype-Resolved 3D Chromatin Maps**

Caspar Gross<sup>1,2,§</sup>, Ramya Potabattula<sup>1,§</sup>, Fubo Cheng<sup>1</sup>, Sarah Leuchtenberg<sup>1</sup>, Hanna Sophie Hartung<sup>1</sup>, Beate Kristmann<sup>1</sup>, Elena Buena-Atienza<sup>1,3</sup>, Nicolas Casadei<sup>1,3</sup>, Stephan Ossowski<sup>1,2,\*</sup> and Olaf Riess<sup>1,3,\*</sup>

<sup>1</sup> Institute of Medical Genetics and Applied Genomics, University of Tübingen, Tübingen, Germany  
<sup>2</sup> Institute for Bioinformatics and Medical Informatics (IBMI), University of Tübingen, Tübingen, Germany  
<sup>3</sup> NGS Competence Center Tübingen (NCCT), Tübingen, Germany  
<sup>§</sup> These authors contributed equally


## Repository structure

```
workflow/       Snakemake workflow (rules, scripts, environments)
data/           Sample and dataset configuration
assembly/       Assembly configuration and QC outputs
doc/            Analysis documentation
```

## Documentation

- [**Datasets**](doc/01_datasets.md): Reference genomes, databases, and published datasets used for benchmarking
- [**Assembly QC**](doc/04_assembly_qc.md): QC metrics and their definitions
- [**Parent-of-Origin**](doc/12_parent_of_origin.md): Methylation-based parent-of-origin assignment

## Quick Start

The analysis is implemented as a Snakemake workflow. If you want to adapt or rerun the analysis, see [workflow/README.md](workflow/README.md) for setup and usage instructions.

### Requirements

- Software: Snakemake ≥ 8.0 with Conda/Mamba
- Compute: GPU resources for HERRO error correction and Dorado basecalling, at least 500GB Ram for Assembly, minimum 64 Threads
- Sequencing data: 3× ultra-long (UL) ONT flowcells and 1× Pore-C flowcell per sample

### Basic usage

```bash
# Configure samples in data/datasets.yml and assembly/assemblies.yml

snakemake --snakefile workflow/Snakefile \
    --use-conda all_samples
```



## Data availability

Generated results, logs, and benchmark files are not stored in this
repository. Together with the phased genome assemblies (FASTA) for all
24 samples, they are archived on Zenodo:
[10.5281/zenodo.23595767](https://doi.org/10.5281/zenodo.23595767).
Raw sequencing reads for this study's own participants are
controlled-access individual-level human genome data, available via
GHGA (study accession GHGAS08644233257767), not on Zenodo.

## Citing this repository

See [CITATION.cff](CITATION.cff). This repository is archived on Zenodo for
citation with a version-specific DOI: **DOI TBD**.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
