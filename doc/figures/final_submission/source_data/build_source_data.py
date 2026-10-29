import pandas as pd
import glob, os, re

ROOT = "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT"
os.chdir(ROOT)
TSV_OUT = "/home/ukt.ad.local/ahgrosc1/.claude/jobs/b85520b6/tmp/source_data_tsv_v2"

sheets = {}

def add(name, df, desc, sources):
    assert name not in sheets, name
    sheets[name] = (df, desc, sources)

def read(path, **kw):
    return pd.read_csv(path, sep="\t" if path.endswith((".tsv", ".txt")) else ",", **kw)

def classify_feature_group(feature_type):
    ft = re.sub(r"_windowed$", "", str(feature_type))
    if re.search(r"Err|MISJOIN", ft, re.IGNORECASE):
        return "Error/Misjoin"
    if re.search(r"COLLAPSE|Col", ft, re.IGNORECASE):
        return "Collapse"
    return ft if ft else "Other"

# ============================= MAIN FIGURES (unchanged from previous build) =============================
add("Fig2A_assembly_metrics", read("doc/tables/qc_samples_table_full.tsv"),
    "Fig 2A: per-sample assembly QC metrics", ["doc/tables/qc_samples_table_full.tsv"])
add("Fig2A_HPRC_comparison", read("doc/tables/qc_hprc_full.tsv"),
    "Fig 2A: HPRC comparison-cohort assembly QC metrics", ["doc/tables/qc_hprc_full.tsv"])

issues_raw = read("doc/tables/for_plots/assembly_issues_all_samples.tsv")
issues_raw["start"] = pd.to_numeric(issues_raw["start"], errors="coerce")
issues_raw["end"] = pd.to_numeric(issues_raw["end"], errors="coerce")
issues_raw = issues_raw.dropna(subset=["start", "end"]).copy()
issues_raw = issues_raw[issues_raw["feature_type"] != "Hap"].copy()
issues_raw["feature_group"] = issues_raw["feature_type"].apply(classify_feature_group)
issues_raw["issue_length"] = issues_raw["end"] - issues_raw["start"]
fig2b = (issues_raw.groupby(["sample", "chr", "feature_group"], as_index=False)
         .agg(n_issues=("issue_length", "size"), affected_length=("issue_length", "sum")))
add("Fig2B_T2T_heatmap_summary", fig2b,
    "Fig 2B / Supp S9: per-sample x chromosome x issue-category summary underlying the ideogram heatmap "
    "(plot further bins into 2Mb windows with per-sample weighting, see 07_figure_1_assembly.ipynb cells 26-27)",
    ["doc/tables/for_plots/assembly_issues_all_samples.tsv (aggregated)"])

issues_unlifted = read("doc/tables/for_plots/assembly_issues_unlifted.tsv")
issues_unlifted = issues_unlifted[(issues_unlifted["feature_type"] != "Hap") &
                                    (issues_unlifted["haplotype"] != "unassigned")].copy()
issues_unlifted["feature_group"] = issues_unlifted["feature_type"].apply(classify_feature_group)
issues_unlifted["issue_length"] = issues_unlifted["end"] - issues_unlifted["start"]
fig2c = (issues_unlifted.groupby(["sample", "feature_type", "feature_group", "source", "haplotype"], as_index=False)
         .agg(n_issues=("issue_length", "size"), affected_length=("issue_length", "sum")))
add("Fig2C_issue_classification", fig2c,
    "Fig 2C: assembly issue counts/length by category, tool and haplotype (dt_issues_classified, "
    "07_figure_1_assembly.ipynb cell 21)",
    ["doc/tables/for_plots/assembly_issues_unlifted.tsv (aggregated)"])

add("Fig3A_ancestry_map", read("doc/tables/samples_metadata.tsv"),
    "Fig 3A: participant birth country / superpopulation", ["doc/tables/samples_metadata.tsv"])
add("Fig3B_global_ancestry", read("doc/tables/for_plots/rfmix_all_chr.msp.tsv"),
    "Fig 3B: per-sample global ancestry proportions (RFMix)", ["doc/tables/for_plots/rfmix_all_chr.msp.tsv"])
add("Fig3C_small_variants_T2T", read("doc/tables/for_plots/small_variant_stats_summary.all_samples.tsv"),
    "Fig 3C: small variant counts, this study's samples", ["doc/tables/for_plots/small_variant_stats_summary.all_samples.tsv"])
add("Fig3C_small_variants_HPRC", read("doc/tables/for_plots/small_variant_stats_summary.hprc.tsv"),
    "Fig 3C: small variant counts, HPRC comparison cohort", ["doc/tables/for_plots/small_variant_stats_summary.hprc.tsv"])
add("Fig3D_SV_stats_T2T", read("doc/tables/for_plots/sv_stats_summary.all_samples.tsv"),
    "Fig 3D: SV counts per chromosome/haplotype, this study's samples", ["doc/tables/for_plots/sv_stats_summary.all_samples.tsv"])
add("Fig3D_SV_stats_HPRC", read("doc/tables/for_plots/sv_stats_summary.hprc.tsv"),
    "Fig 3D: SV counts per chromosome/haplotype, HPRC comparison cohort", ["doc/tables/for_plots/sv_stats_summary.hprc.tsv"])
add("Fig4C_contact_counts", read("analysis_other/dip3d/combined_dip3d_stats.tsv"),
    "Fig 4C: Pore-C contact counts per sample", ["analysis_other/dip3d/combined_dip3d_stats.tsv"])

cen = read("doc/tables/for_plots/cenmap_single_complete_correct_cens.all_samples.tsv")
add("Fig5A_centromeres_per_chr", cen, "Fig 5A: centromere calls, by chromosome",
    ["doc/tables/for_plots/cenmap_single_complete_correct_cens.all_samples.tsv"])
add("Fig5B_centromeres_per_sample", cen.copy(), "Fig 5B: centromere calls, by sample",
    ["doc/tables/for_plots/cenmap_single_complete_correct_cens.all_samples.tsv"])

# ============================= SUPPLEMENTARY FIGURES (corrected S-numbering) =============================

# S1 - sequencing yield, all flow cells
add("S1_flowcell_yield", read("doc/tables/flowcell_stats.tsv"),
    "Supp S1: per-flowcell sequencing yield/N50 across all datasets",
    ["doc/tables/flowcell_stats.tsv"])

# S2 - UL read length optimization over time. Rebuilt to match the original figure (extracted from the
# supplement docx by the user): scatter of Read N50 (kbp) vs Flowcell yield (Gbp), points colored by
# sequencing date on a rainbow/turbo scale. Not found in any commit of the dev repo's git history
# (never version-controlled); dates parsed from raw MinKNOW run-folder timestamps. T2T_HG002 flowcells
# (dated Feb 2026) excluded -- they postdate the original figure, whose color scale tops out ~Sep 2025.
add("S2_ul_read_length_by_date", read("doc/tables/for_plots/ul_read_length_by_date.tsv"),
    "Supp S2: UL flow cell read N50 vs flowcell yield, colored by sequencing date (dates parsed from "
    "raw MinKNOW run-folder names, e.g. data/raw/.../<flowcell_id>/<YYYYMMDD_HHMM>_.../; T2T_HG002 "
    "flowcells excluded as out of scope for this timeline)",
    ["doc/tables/for_plots/ul_read_length_by_date.tsv"])

# S3 - assembly input read statistics
# NOTE: uses the non-"_all" table. The "_all" table (used in the first pass) has since been updated
# with post-publication samples (T2T20-T2T53) that are not part of this study; the plain table matches
# the publication cohort only, same as the figure already published in the public T2T_ONT repo
# (doc/figures/input_qc_read_stats.png).
add("S3_assembly_input_read_stats", read("doc/tables/for_plots/assembly_input_read_stats.tsv"),
    "Supp S3: HQ/UL input read statistics per sample, publication cohort only",
    ["doc/tables/for_plots/assembly_input_read_stats.tsv"])

# S4 - runtime -> pointer only (multi-hundred log aggregation, no single saved table)

# S5 - detailed chromosome assembly classification (per-contig)
add("S5_t2t_grid_detailed", read("doc/tables/for_plots/T2T_scaffolds.both.all.txt"),
    "Supp S5: per-contig T2T/scaffold classification with gap statistics, all samples/haplotypes",
    ["doc/tables/for_plots/T2T_scaffolds.both.all.txt"])

# S6 - assembled chrY length (simplified: contig/scaffold classification only, chrY rows)
t2t_scaffolds = read("doc/tables/for_plots/T2T_scaffolds.both.all.txt")
s6 = t2t_scaffolds[t2t_scaffolds["chromosome"] == "chrY"].copy()
add("S6_chrY_length", s6,
    "Supp S6: assembled chrY length and T2T/scaffold classification per sample (simplified -- the "
    "published figure additionally overlays gap/satellite tracks per sample, sourced from "
    "assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed and "
    "analysis_other/repeatmasker/{sample}/rm_summary/{sample}_filtered_satellites.bed, not flattened here; "
    "see workflow/scripts/59_plot_Y_karyogram.R header for the full input list)",
    ["doc/tables/for_plots/T2T_scaffolds.both.all.txt (filtered chromosome == chrY)"])

# S7 - extended QC metrics vs HPRC: identical data to Fig2A sheets, referenced not duplicated

# S8 - Flagger issue length by sample (filtered view of Fig2C source)
s8 = fig2c[fig2c["source"] == "Flagger"].copy()
add("S8_flagger_issue_length", s8,
    "Supp S8: Flagger-detected assembly issue counts/length by sample and category (same source as "
    "Fig2C_issue_classification, filtered to source == Flagger)",
    ["doc/tables/for_plots/assembly_issues_unlifted.tsv (aggregated, Flagger only)"])

# S9 - issues heatmap WITHOUT the "Known Difficult Regions" panel (the published figure does not
# include it -- that panel relies on external literature-sourced coordinates the user does not want
# referenced/included). Reuses the same data as Fig2B_T2T_heatmap_summary; no separate sheet needed.

# S10 - multicopy gene distance (gaps + contig ends)
gap_files = glob.glob("analysis_other/multicopy_genes/*/gap_to_multicopy_gene_distance.tsv")
gap_all = []
for f in gap_files:
    sample = os.path.basename(os.path.dirname(f))
    d = pd.read_csv(f, sep="\t")
    d.insert(0, "sample", sample)
    gap_all.append(d)
gap_df = pd.concat(gap_all, ignore_index=True)
add("S10_gap_to_gene_distance", gap_df,
    "Supp S10: distance from assembly gaps to nearest multi-copy gene locus, all samples",
    [f"analysis_other/multicopy_genes/<sample>/gap_to_multicopy_gene_distance.tsv ({len(gap_files)} samples)"])

# S10 contig-end panel dropped from the workbook: at ~189k rows its per-cell XML bloats the sheet
# to ~123MB uncompressed, which was making Excel flag the whole file as needing repair. Pointer only;
# raw per-sample files are archived on Zenodo.
ce_files = glob.glob("analysis_other/multicopy_genes/*/contig_ends_to_multicopy_gene_distance.tsv")

# S11 - small contig satellite classification
add("S11_small_contig_sat_class", read("doc/tables/for_plots/small_contig_satellite_class.all_samples.tsv"),
    "Supp S11: small-contig satellite classification, all samples",
    ["doc/tables/for_plots/small_contig_satellite_class.all_samples.tsv"])

# S12 - QV across k-mers (A) and strategies (B)
qv = read("doc/tables/for_plots/all_merqury_qv.tsv")
s12a = qv[(qv["modifier"].isna() | (qv["modifier"] == "")) & (qv["haplotype"] != "diploid")].copy()
add("S12A_qv_by_kmer", s12a,
    "Supp S12A: consensus QV by k-mer size (21/31/37), default assemblies, per haplotype",
    ["doc/tables/for_plots/all_merqury_qv.tsv (filtered modifier == default, haplotype != diploid)"])
s12b = qv[(qv["sample"].isin(["T2T00", "T2T04", "T2T03"])) & (qv["haplotype"] == "diploid")].copy()
add("S12B_qv_by_strategy", s12b,
    "Supp S12B: consensus QV comparison across strategies (default Pore-C / trio / APK-polished) for "
    "T2T00, T2T04, T2T03, diploid",
    ["doc/tables/for_plots/all_merqury_qv.tsv (filtered sample in T2T00/T2T04/T2T03, haplotype == diploid)"])

# S13 - HG002 error classification + QV by variant type (same source table, both panels)
add("S13_HG002_error_class_and_qv", read("analysis_other/GQC/T2T_HG002/gqc_error_classification/summary.tsv"),
    "Supp S13: HG002 assembly error counts by genome-context category (panel A) and QV by variant "
    "type for reliable categories (panel B) -- same source table for both panels "
    "(16_additional_qc.ipynb cells 26-27)",
    ["analysis_other/GQC/T2T_HG002/gqc_error_classification/summary.tsv"])

# S14 - HG002 GQC native output (generalstats + homopolymer length summary as flat-table proxies;
# the actual A-D panels are GQC's own rendered PDF plots)
add("S14_HG002_GQC_generalstats",
    pd.read_csv("analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.generalstats.txt", sep="\t", header=None, names=["metric", "value"]),
    "Supp S14: GQC assemblybench summary statistics for HG002 (direct GQC tool output; panels A-D of "
    "the published figure are GQC's own rendered plots -- "
    "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.{continuitystats,mononuc_accuracy,"
    "discrepancystats,indelsizestats}.HG002v1.1.pdf)",
    ["analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.generalstats.txt"])
add("S14_HG002_hp_length_summary", read("analysis_other/GQC/T2T_HG002/gqc_error_classification/hp_length_summary.tsv"),
    "Supp S14: homopolymer-length error summary for HG002 (direct GQC tool output)",
    ["analysis_other/GQC/T2T_HG002/gqc_error_classification/hp_length_summary.tsv"])

# S15 dropped from the workbook: filtered to sample T2T04 this is still ~441k rows, bloating its sheet
# to ~300MB uncompressed XML -- the main cause of Excel's repair prompt. Pointer only; raw phase block
# data archived on Zenodo.

# S16 - assembler/phasing QC method comparison (12 metrics incl. genome_completeness, MMC, MSC).
# Identified by the user directly from 06_assembly_types.ipynb cell 5 (plot_qc_method_comparison /
# trioeval_qc_comparison.png). Reproduced here by re-running that cell's exact R processing
# (process_qc_table(selection = metrics_phased) on assembly/qc/qc_trio.tsv) via the r-notebooks conda
# env, since this long-format table isn't saved anywhere by the notebook itself; the result is saved to
# doc/tables/for_plots/qc_trio_phased_metrics.tsv for reuse.
add("S16_qc_method_comparison", read("doc/tables/for_plots/qc_trio_phased_metrics.tsv"),
    "Supp S16: assembly/phasing QC metrics (length, T2T contigs/scaffolds, completeness, MMC, MSC, "
    "gaps, misjoins) by assembler (Verkko/Hifiasm) x phasing method (Trio/Pore-C), samples T2T00/T2T04",
    ["doc/tables/for_plots/qc_trio_phased_metrics.tsv (reproduced from assembly/qc/qc_trio.tsv via "
     "06_assembly_types.ipynb cells 3-5 processing)"])

# S17 - switch/hamming error rates
add("S17_switch_hamming_errors", read("analysis_other/trioeval/trioeval_summary.tsv"),
    "Supp S17: switch and Hamming error rates by phasing method, sample and haplotype",
    ["analysis_other/trioeval/trioeval_summary.tsv"])

# S18 - parent-of-origin assignment across cohort (concatenated per-sample scores)
pofo_samples = ["T2T00", "T2T01", "T2T02", "T2T03", "T2T04", "T2T05", "T2T06", "T2T07", "T2T08", "T2T09",
                 "T2T10", "T2T11", "T2T12", "T2T13", "T2T15", "T2T16", "T2T17", "T2T18", "T2T19", "T2T_HG002"]
pofo_all = []
for s in pofo_samples:
    f = f"analysis_other/parent_of_origin/{s}/{s}_PofO_Scores.tsv"
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t")
        d.insert(0, "sample", s)
        pofo_all.append(d)
pofo_df = pd.concat(pofo_all, ignore_index=True)
add("S18_pofo_cohort", pofo_df,
    "Supp S18: per-chromosome parent-of-origin assignment scores across the cohort "
    f"({len(pofo_all)} of {len(pofo_samples)} samples found)",
    ["analysis_other/parent_of_origin/<sample>/<sample>_PofO_Scores.tsv (concatenated)"])

# S19 - Ti/Tv ratio (A) + heterozygous proportion (B)
small_var_t2t = read("doc/tables/for_plots/small_variant_stats_summary.all_samples.tsv")
small_var_hprc = read("doc/tables/for_plots/small_variant_stats_summary.hprc.tsv")
small_var_t2t["dataset"] = "Samples"
small_var_hprc["dataset"] = "HPRC"
small_var_combined = pd.concat([small_var_t2t, small_var_hprc], ignore_index=True)
s19a = small_var_combined[(small_var_combined["category"] == "ti_tv") &
                           (small_var_combined["subcategory"] == "transitions")].copy()
add("S19A_titv_ratio", s19a,
    "Supp S19A: transition/transversion (Ti/Tv) ratio per sample (percentage column stores the ratio)",
    ["doc/tables/for_plots/small_variant_stats_summary.{all_samples,hprc}.tsv (filtered category == ti_tv)"])
s19b = small_var_combined[(small_var_combined["category"] == "zygosity") &
                           (small_var_combined["var_type"].isin(["SNV", "DEL", "INS"])) &
                           (small_var_combined["subcategory"] == "HET")].copy()
add("S19B_het_proportion", s19b,
    "Supp S19B: proportion of heterozygous variants per genome, by variant type",
    ["doc/tables/for_plots/small_variant_stats_summary.{all_samples,hprc}.tsv (filtered category == zygosity, subcategory == HET)"])

# S20 - chr9 SV density (CORRECTED: binned mean-per-sample, matching the actual plotted data,
# not the raw per-variant detail table used in the previous build)
CENSAT_START, CENSAT_END, BIN_SIZE = 40.0e6, 81.7e6, 2.0e6
sv_details_t2t = read("doc/tables/for_plots/sv_stats_details.all_samples.tsv")
chr9 = sv_details_t2t[sv_details_t2t["chrom"] == "chr9"].copy()
chr9["bin_start"] = (chr9["pos"] // BIN_SIZE) * BIN_SIZE
chr9["region"] = chr9["pos"].apply(
    lambda p: "Pericentromeric heterochromatin (9q12)" if CENSAT_START <= p <= CENSAT_END
    else ("9p arm" if p < CENSAT_START else "9q arm"))
counts = chr9.groupby(["sample", "bin_start", "region"]).size().reset_index(name="n")
chr9_density = (counts.groupby(["bin_start", "region"], as_index=False)
                 .agg(mean_count=("n", "mean")))
chr9_density["bin_mid_Mb"] = (chr9_density["bin_start"] + BIN_SIZE / 2) / 1e6
add("S20_chr9_sv_density", chr9_density,
    "Supp S20: mean SV count per 2Mb bin along chr9, by region, averaged across samples "
    "(exact aggregation plotted in 08_figure_2_variations.ipynb cell 17; corrects the previous build, "
    "which used the raw per-variant detail table instead of this binned mean)",
    ["doc/tables/for_plots/sv_stats_details.all_samples.tsv (filtered chrom == chr9, binned 2Mb, mean per sample)"])

# S21 - PCA of study cohort
pca = pd.read_csv("analysis_other/ancestry/pca/merge_ref_samples.eigenvec", sep="\t", comment="#", header=None)
pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca.shape[1] - 1)]
psam = pd.read_csv("analysis_other/ancestry/plink/merge_ref_samples.psam", sep="\t", header=0)
psam.columns = [c.lstrip("#") for c in psam.columns]
psam = psam[["IID", "SUPERPOP"]]
pca_merged = pca.merge(psam, on="IID", how="left")
add("S21_PCA_cohort", pca_merged,
    "Supp S21: principal components (PLINK2 eigenvec) of phased small variants merged with the "
    "T2T-lifted 1000 Genomes reference panel, with superpopulation labels",
    ["analysis_other/ancestry/pca/merge_ref_samples.eigenvec", "analysis_other/ancestry/plink/merge_ref_samples.psam"])

# S22 - Pore-C sequencing statistics (filtered view of the same table as S3; same "_all" vs
# publication-cohort fix applied here as for S3)
s22 = read("doc/tables/for_plots/assembly_input_read_stats.tsv")
if "type" in s22.columns:
    s22 = s22[s22["type"] == "POREC"].copy()
add("S22_porec_input_stats", s22,
    "Supp S22: Pore-C input read statistics per sample, publication cohort only (same source table as "
    "S3, filtered to Pore-C)",
    ["doc/tables/for_plots/assembly_input_read_stats.tsv (filtered type == POREC)"])

# S23 - per-flowcell Pore-C QC (wf-pore-c / MultiQC)
add("S23_porec_flowcell_qc", read("analysis_other/porec/PoreC-QC_multiqc_report_data/multiqc_general_stats.txt"),
    "Supp S23: wf-pore-c quality metrics per flowcell/haplotype (MultiQC general stats)",
    ["analysis_other/porec/PoreC-QC_multiqc_report_data/multiqc_general_stats.txt"])

# ============================= Pointer-only panels =============================
pointer_rows = [
    {"figure": "Fig 1 (overview)", "reason": "Schematic/infographic, no plotted dataset", "where": "n/a (schematic)"},
    {"figure": "Fig 4A/B (Hi-C contact matrices)", "reason": "Multi-track genome-browser visualization",
     "where": "Raw Pore-C contact data on Zenodo (10.5281/zenodo.23595767); code in 09_figure_3_3D.ipynb, 09_figure_3_fanc.ipynb"},
    {"figure": "Fig 5C/D (HOR structure / self-identity dot plots)", "reason": "BED/matrix-based genome tracks, custom loader",
     "where": "Raw centromere annotation on Zenodo (10.5281/zenodo.23595767); code in 11_figure_4_cenplot.ipynb"},
    {"figure": "Supp S4 (computational resource requirements / runtime)",
     "reason": "Aggregated over several hundred per-run log/benchmark files, no single saved source table",
     "where": "Raw logs archived on Zenodo (10.5281/zenodo.23595767), directory runtimes/; code in 17_runtimes.ipynb"},
    {"figure": "Supp S10 (contig-end to multicopy-gene distance, panel B)",
     "reason": f"~{len(ce_files)}-sample per-event table, ~189k rows -- too large for an Excel sheet "
               "(uncompressed worksheet XML ~123MB, was causing the file to need repair)",
     "where": "Raw per-sample files (analysis_other/multicopy_genes/<sample>/contig_ends_to_multicopy_gene_distance.tsv) "
              "archived on Zenodo (10.5281/zenodo.23595767); code in 16_additional_qc.ipynb"},
    {"figure": "Supp S15 (phase block karyogram, T2T04, 5 methods)",
     "reason": "~441k rows even filtered to one sample -- too large for an Excel sheet "
               "(uncompressed worksheet XML ~300MB, the main cause of the file needing repair)",
     "where": "Raw phase block calls (analysis_other/trioeval/trioeval_phaseblocks.tsv) archived on "
              "Zenodo (10.5281/zenodo.23595767); code in 06_assembly_types.ipynb"},
]
pointer_df = pd.DataFrame(pointer_rows)

# ============================= Write outputs =============================
index_rows = []
for name, (df, desc, sources) in sheets.items():
    index_rows.append({"sheet": name, "description": desc, "rows": len(df), "columns": len(df.columns),
                        "source_file(s)_in_dev_repo": "; ".join(sources)})
index_df = pd.DataFrame(index_rows)

print(f"{len(sheets)} sheets built, writing...")
os.makedirs(TSV_OUT, exist_ok=True)
for name, (df, desc, sources) in sheets.items():
    df.to_csv(f"{TSV_OUT}/{name}.tsv", sep="\t", index=False)
index_df.to_csv(f"{TSV_OUT}/00_README_index.tsv", sep="\t", index=False)
pointer_df.to_csv(f"{TSV_OUT}/00_README_pointers_not_included.tsv", sep="\t", index=False)

with pd.ExcelWriter("doc/publication/revision2/source_data/SourceData.xlsx", engine="xlsxwriter") as writer:
    index_df.to_excel(writer, sheet_name="README_index", index=False)
    pointer_df.to_excel(writer, sheet_name="README_not_included", index=False)
    for name, (df, desc, sources) in sheets.items():
        df.to_excel(writer, sheet_name=name[:31], index=False)

print("DONE")
print(index_df[["sheet", "rows", "columns"]].to_string())
