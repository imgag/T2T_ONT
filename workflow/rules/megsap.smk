"""
Snakemake workflow for megSAP longread analysis.
Steps:
  1. Copy input BAM files to megSAP sample folders (UL runs only)
  2. Merge multiple UL runs per biological sample
  3. Run megSAP analyze_longread.php on the final merged sample
"""

import os
from glob import glob

##############################################################################
# Configuration
##############################################################################

MEGSAP_DIR     = "megSAP"
BASECALLED_DIR = "data/basecalled/sup"
MEGSAP_PHP     = "/mnt/storage2/users/ahgrosc1/dev/megsap/src/Pipelines/analyze_longread.php"
MERGE_PHP      = "/mnt/storage2/megSAP/pipeline/src/Tools/merge_samples.php"
THREADS        = 32
ASSOC_TABLE    = "doc/tables/flowcell_biological_sample.tsv"

# Only this processing system is used for megSAP analysis
UL_PROCESSING_SYSTEM = "UL-LR-ONT-SQK-ULK114"

# ---- Provide the list of biological sample names (name_external) to process ----
TARGET_SAMPLES = [
    "T2T20",
    "T2T21"
]

##############################################################################
# Helper: parse association table — UL runs only
# TSV columns: name_ngsd  name_external  project_name  run_flowcell_id  run_id  processing_system_name
##############################################################################

def parse_ul_runs(assoc_table_path, target_samples, ul_ps=UL_PROCESSING_SYSTEM):
    """
    Parse the association table and return only UL-LR-ONT-SQK-ULK114 runs
    for each biological sample in target_samples.

    Returns:
        ul_runs: dict  biological_sample -> [(ngsd_id, run_id), ...]  (sorted by ngsd_id)
    """
    ul_runs = {}

    with open(assoc_table_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("name_ngsd"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue

            name_ngsd, name_external, _, _, run_id, processing_system = parts[:6]

            # Skip anything that is not a UL run
            if processing_system != ul_ps:
                continue

            # Skip samples not in the target list
            if name_external not in target_samples:
                continue

            # Strip leading '#' from run_id (e.g. '#06113' -> '06113')
            run_id_clean = run_id.lstrip("#").strip()

            ul_runs.setdefault(name_external, []).append((name_ngsd, run_id_clean))

    # Sort by run_id (numeric) for deterministic merge order
    for k in ul_runs:
        ul_runs[k] = sorted(ul_runs[k], key=lambda x: int(x[1]))

    return ul_runs


UL_RUNS = parse_ul_runs(ASSOC_TABLE, TARGET_SAMPLES)

# Validate that every target sample has at least one UL run
for s in TARGET_SAMPLES:
    if s not in UL_RUNS:
        raise ValueError(
            f"No UL-LR-ONT-SQK-ULK114 runs found for target sample '{s}' "
            f"in {ASSOC_TABLE}. Check TARGET_SAMPLES or the association table."
        )

##############################################################################
# Helper: locate source BAM folder using run_id as unique anchor
##############################################################################

def get_source_bam(ngsd_id, run_id):
    """
    Find the BAM file for a given ngsd sample ID under BASECALLED_DIR.
    Matches the folder solely by run_id (e.g. '06113'), which is unique
    and handles all folder naming variants:
      25006LRa212_T2T21_06113
      25006LRa212_T2T21_UL_06113
    """
    # Match any subfolder ending with the run_id
    folder_pattern = os.path.join(BASECALLED_DIR, f"*{run_id}")
    folders = glob(folder_pattern)

    if not folders:
        raise ValueError(
            f"No basecalled folder found for {ngsd_id} (run_id={run_id}) "
            f"matching pattern: {folder_pattern}"
        )
    if len(folders) > 1:
        raise ValueError(
            f"Multiple basecalled folders matched run_id={run_id}: "
            f"{folders} — run_id is not unique."
        )

    folder = folders[0]
    # Use the base ngsd_id (before underscore) to match the BAM file
    ngsd_base = ngsd_id.split("_")[0]
    bam_hits = glob(os.path.join(folder, f"{ngsd_base}_*.bam"))

    if not bam_hits:
        raise ValueError(
            f"No BAM file found in folder {folder} matching {ngsd_base}_*.bam"
        )

    return bam_hits[0]




def dest_bam_path(ngsd_id):
    """Destination BAM path inside the megSAP sample folder."""
    return os.path.join(MEGSAP_DIR, f"Sample_{ngsd_id}", f"{ngsd_id}.mod.unmapped.bam")


# Build flat list of copy jobs: (ngsd_id, run_id, src, dst)
# Only UL runs are included — PoreC / APK / LR runs are excluded by parse_ul_runs()
copy_jobs = []
for bio_sample, id_pairs in UL_RUNS.items():
    for ngsd_id, run_id in id_pairs:
        src = get_source_bam(ngsd_id, run_id)
        dst = dest_bam_path(ngsd_id)
        copy_jobs.append((ngsd_id, run_id, src, dst))

# Map dst -> src for use in rule input lambdas
dst_to_src = {dst: src for (_, _, src, dst) in copy_jobs}

# Final merged sample per biological sample = last ngsd_id after sorting by run_id
final_ul_id = {bio: pairs[-1][0] for bio, pairs in UL_RUNS.items()}


##############################################################################
# Target rule
##############################################################################

rule all:
    input:
        expand(
            os.path.join(MEGSAP_DIR, "{bio_sample}_analysis.done"),
            bio_sample=TARGET_SAMPLES
        )

##############################################################################
# Step 1 – Copy UL BAM files into megSAP sample folders
##############################################################################

rule copy_bam:
    """Copy a single UL flowcell BAM into the megSAP sample folder."""
    input:
        lambda wc: dst_to_src[
            dest_bam_path(wc.ngsd_id)
        ]
    output:
        os.path.join(
            MEGSAP_DIR,
            "Sample_{ngsd_id}",
            "{ngsd_id}.mod.unmapped.bam"
        )
    log:
        "logs/copy_bam/{ngsd_id}.log"
    shell:
        """
        mkdir -p $(dirname {output})
        cp {input} {output} 2> {log}
        echo "Copied {input} -> {output}" >> {log}
        """

##############################################################################
# Step 2 – Merge multiple UL runs of the same biological sample
##############################################################################

def get_all_ul_bams(bio_sample):
    """Return all destination BAMs for the UL runs of a biological sample."""
    return [
        dest_bam_path(ngsd_id)
        for ngsd_id, _ in UL_RUNS[bio_sample]
    ]

rule merge_runs:
    """
    Sequentially merge all UL runs into the last sample using merge_samples.php.
    No-op when only a single run exists.
    """
    input:
        lambda wc: get_all_ul_bams(wc.bio_sample)
    output:
        touch(os.path.join(MEGSAP_DIR, "{bio_sample}_merged.done"))
    params:
        ids        = lambda wc: [ngsd_id for ngsd_id, _ in UL_RUNS[wc.bio_sample]],
        merge_php  = MERGE_PHP,
        megsap_dir = MEGSAP_DIR,
    log:
        "logs/merge/{bio_sample}.log"
    shell:
        """
        set -euo pipefail
        ids=({params.ids})
        if [[ ${{#ids[@]}} -le 1 ]]; then
            echo "Single run — no merging needed." | tee {log}
            exit 0
        fi
        pushd {params.megsap_dir} > /dev/null
        for ((i=0; i < ${{#ids[@]}}-1; i++)); do
            current="${{ids[$i]}}"
            next="${{ids[$i+1]}}"
            echo "Merging ${{current}} into ${{next}}" | tee -a ../{log}
            php {params.merge_php} -ps "${{current}}" -into "${{next}}" >> ../{log} 2>&1
        done
        popd > /dev/null
        """

##############################################################################
# Step 3 – Run megSAP analyze_longread.php on the final merged UL sample
##############################################################################

rule analyze_longread:
    """Run megSAP longread analysis on the final (merged) UL sample."""
    input:
        os.path.join(MEGSAP_DIR, "{bio_sample}_merged.done")
    output:
        touch(os.path.join(MEGSAP_DIR, "{bio_sample}_analysis.done"))
    params:
        final_id   = lambda wc: final_ul_id[wc.bio_sample],
        megsap_php = MEGSAP_PHP,
        megsap_dir = MEGSAP_DIR,
        threads    = THREADS,
    log:
        "logs/analyze/{bio_sample}.log"
    shell:
        """
        set -euo pipefail
        cd {params.megsap_dir}
        php {params.megsap_php} \
            -folder Sample_{params.final_id} \
            -name   {params.final_id} \
            -threads {params.threads} \
            > ../{log} 2>&1
        """