#!/usr/bin/env python3
"""
58_aggregate_porec_stats.py

Aggregate pairtools stats files across all samples and haplotypes into a
single summary TSV.  For each input stats file, the sample ID, haplotype
and pipeline of origin are inferred from the directory name:

  analysis_other/porec/{sample}/       -> haplotype=unphased  pipeline=wf-pore-c
  analysis_other/porec/{sample}.hp1/   -> haplotype=hp1       pipeline=dip3d
  analysis_other/porec/{sample}.hp2/   -> haplotype=hp2       pipeline=dip3d

Metrics extracted (all raw counts + derived rates):
  total, total_unmapped, total_single_sided_mapped, total_mapped,
  total_dups, total_nodups, cis, trans,
  cis_1kb+, cis_2kb+, cis_10kb+, cis_40kb+,
  frac_cis, frac_cis_1kb+, frac_cis_2kb+, frac_cis_10kb+, frac_cis_40kb+,
  frac_dups, complexity_naive,
  mapping_rate, single_sided_rate, unmapped_rate,
  convergence_dist

Usage
-----
  python 58_aggregate_porec_stats.py \
      --stats file1.pairs.stats.txt [file2 ...] \
      --output doc/tables/for_plots/porec_stats_summary.tsv
"""

import argparse
import os
import re
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Metric keys extracted verbatim from the stats file (key → output column)
# ---------------------------------------------------------------------------
RAW_METRICS = {
    "total":                              "total",
    "total_unmapped":                     "total_unmapped",
    "total_single_sided_mapped":          "total_single_sided_mapped",
    "total_mapped":                       "total_mapped",
    "total_dups":                         "total_dups",
    "total_nodups":                       "total_nodups",
    "cis":                                "cis",
    "trans":                              "trans",
    "cis_1kb+":                           "cis_1kb",
    "cis_2kb+":                           "cis_2kb",
    "cis_10kb+":                          "cis_10kb",
    "cis_40kb+":                          "cis_40kb",
    "summary/frac_cis":                   "frac_cis",
    "summary/frac_cis_1kb+":             "frac_cis_1kb",
    "summary/frac_cis_2kb+":             "frac_cis_2kb",
    "summary/frac_cis_10kb+":            "frac_cis_10kb",
    "summary/frac_cis_40kb+":            "frac_cis_40kb",
    "summary/frac_dups":                  "frac_dups",
    "summary/complexity_naive":           "complexity_naive",
    "summary/dist_freq_convergence/convergence_dist": "convergence_dist",
}


def parse_identity(stats_path: str) -> tuple[str, str, str]:
    """
    Derive (sample, haplotype, pipeline) from a stats file path.

    Expected path structure (relative or absolute):
      .../porec/{dataset}/pairs/{dataset}.pairs.stats.txt

    Where {dataset} is one of:
      {sample}        → unphased, wf-pore-c
      {sample}.hp1    → hp1,      dip3d
      {sample}.hp2    → hp2,      dip3d
    """
    # Use the parent-of-parent directory as the dataset identifier
    dataset = os.path.basename(os.path.dirname(os.path.dirname(stats_path)))

    m = re.fullmatch(r"(.+)\.(hp[12])", dataset)
    if m:
        sample, haplotype = m.group(1), m.group(2)
        pipeline = "dip3d"
    else:
        sample = dataset
        haplotype = "unphased"
        pipeline = "wf-pore-c"

    return sample, haplotype, pipeline


def parse_stats_file(path: str) -> dict:
    """
    Parse a pairtools .pairs.stats.txt file and return a flat key→value dict.
    Values are cast to int or float where possible.
    """
    data = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key, value = parts[0], parts[1]
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
            data[key] = value
    return data


def build_record(stats_path: str) -> dict | None:
    """
    Build one output row for a single stats file.
    Returns None and logs a warning if the file does not exist.
    """
    if not os.path.isfile(stats_path):
        print(f"[WARNING] file not found, skipping: {stats_path}", file=sys.stderr)
        return None

    sample, haplotype, pipeline = parse_identity(stats_path)
    raw = parse_stats_file(stats_path)

    record = {
        "sample":    sample,
        "haplotype": haplotype,
        "pipeline":  pipeline,
    }

    # Extract defined metrics
    for key, col in RAW_METRICS.items():
        val = raw.get(key, None)
        # pairtools reports "inf" as a string for complexity_naive when no dups
        if isinstance(val, str) and val.lower() in ("inf", "nan"):
            val = float(val.lower().replace("inf", "inf"))
        record[col] = val

    # Derived rates (guard against zero denominator)
    total = raw.get("total", 0) or 0
    if total > 0:
        record["mapping_rate"]      = raw.get("total_mapped", 0) / total
        record["single_sided_rate"] = raw.get("total_single_sided_mapped", 0) / total
        record["unmapped_rate"]     = raw.get("total_unmapped", 0) / total
    else:
        record["mapping_rate"] = record["single_sided_rate"] = record["unmapped_rate"] = None

    return record


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate pairtools stats files into a single TSV table."
    )
    parser.add_argument(
        "--stats", nargs="+", required=True, metavar="FILE",
        help="One or more .pairs.stats.txt files to aggregate."
    )
    parser.add_argument(
        "--output", required=True, metavar="TSV",
        help="Output path for the aggregated TSV table."
    )
    args = parser.parse_args()

    records = []
    for path in args.stats:
        rec = build_record(path)
        if rec is not None:
            records.append(rec)

    if not records:
        print("[ERROR] No valid stats files found.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(records)

    # Sort by sample → haplotype (unphased first, then hp1, hp2)
    hp_order = {"unphased": 0, "hp1": 1, "hp2": 2}
    df["_hp_order"] = df["haplotype"].map(hp_order).fillna(99)
    df = df.sort_values(["sample", "_hp_order"]).drop(columns=["_hp_order"])
    df = df.reset_index(drop=True)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Written {len(df)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
