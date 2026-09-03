"""
GQC assembly error classification using GIAB HG002v1.1 stratifications.

Annotation layers — all in HG002v1.1 chr_MATERNAL / chr_PATERNAL coordinates:
  cds              JHU Liftoff v0.6 RefSeq CDS lifted onto HG002v1.1
  satellite        GIAB HG002 satellites (cenSat-derived, with 5bp slop)
  segdup           GIAB HG002 segmental duplications >10 kb
  tandem_repeat    GIAB HG002 AllTandemRepeats (excludes homopolymers)
  homopolymer      GIAB HG002 AllHomopolymers ≥7bp + imperfect ≥11bp
  mhc              GIAB HG002 MHC region
  kir              GIAB HG002 KIR region
  vdj              GIAB HG002 VDJ loci
"""

# Map each GIAB category to its (subdirectory, filename suffix after HG002.{hap}_)
_GIAB_STRATS = {
    "tandem_repeat": ("LowComplexity",         "AllTandemRepeats.bed.gz"),
    "homopolymer":   ("LowComplexity",         "AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz"),
    "satellite":     ("LowComplexity",         "satellites_slop5.bed.gz"),
    "segdup":        ("SegmentalDuplications", "segdups_gt10kb.bed.gz"),
    "mhc":           ("OtherDifficult",        "MHC.bed.gz"),
    "kir":           ("OtherDifficult",        "KIR.bed.gz"),
    "vdj":           ("OtherDifficult",        "VDJ.bed.gz"),
}

# Map each intersect category to its annotation BED
_INTERSECT_ANNOTS = {
    "cds":           "data/ref/hg002v1.1.cds_jhu_liftoff.bed",
    "tandem_repeat": config["gqc_error_classification"]["giab_strat_merged"] + "/tandem_repeat.bed",
    "homopolymer":   config["gqc_error_classification"]["giab_strat_merged"] + "/homopolymer.bed",
    "satellite":     config["gqc_error_classification"]["giab_strat_merged"] + "/satellite.bed",
    "segdup":        config["gqc_error_classification"]["giab_strat_merged"] + "/segdup.bed",
    "mhc":           config["gqc_error_classification"]["giab_strat_merged"] + "/mhc.bed",
    "kir":           config["gqc_error_classification"]["giab_strat_merged"] + "/kir.bed",
    "vdj":           config["gqc_error_classification"]["giab_strat_merged"] + "/vdj.bed",
}


# ── Target ────────────────────────────────────────────────────────────────────

rule all_gqc_error_classification:
    input:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/summary.tsv",
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/hp_length_summary.tsv",


# ── CDS annotation BED ───────────────────────────────────────────────────────

rule hg002_cds_bed:
    """Extract CDS features from the JHU Liftoff v0.6 GFF3 files into a merged BED."""
    input:
        mat_gff = config["gqc_error_classification"]["cds_mat_gff"],
        pat_gff = config["gqc_error_classification"]["cds_pat_gff"],
    output:
        "data/ref/hg002v1.1.cds_jhu_liftoff.bed",
    log:
        "logs/gqc/hg002_cds_bed.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        {{ zcat {input.mat_gff}; zcat {input.pat_gff}; }} \
            | grep -v '^#' \
            | awk '$3=="CDS" {{print $1"\t"$4-1"\t"$5}}' \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > {output} 2>{log}
        """


# ── GIAB stratification BEDs: merge mat + pat per category ───────────────────

rule giab_merge_haplotypes:
    """Merge maternal and paternal GIAB stratification BEDs into one merged BED."""
    input:
        mat = lambda wc: (
            config["gqc_error_classification"]["giab_strat_mat"]
            + f"/{_GIAB_STRATS[wc.strat][0]}/HG002.mat_{_GIAB_STRATS[wc.strat][1]}"
        ),
        pat = lambda wc: (
            config["gqc_error_classification"]["giab_strat_pat"]
            + f"/{_GIAB_STRATS[wc.strat][0]}/HG002.pat_{_GIAB_STRATS[wc.strat][1]}"
        ),
    output:
        config["gqc_error_classification"]["giab_strat_merged"] + "/{strat}.bed",
    wildcard_constraints:
        strat = "|".join(_GIAB_STRATS),
    log:
        "logs/gqc/giab_merge_{strat}.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        {{ zcat {input.mat}; zcat {input.pat}; }} \
            | cut -f1-3 \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > {output} 2>{log}
        """


# ── Per-category error intersections ─────────────────────────────────────────

rule gqc_intersect:
    """Intersect GQC assembly errors against one annotation BED; output overlapping errors (-u)."""
    input:
        errors     = "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.errortype.HG002v1.1.bed",
        annotation = lambda wc: _INTERSECT_ANNOTS[wc.category],
    output:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/{category}.bed",
    wildcard_constraints:
        category = "|".join(_INTERSECT_ANNOTS),
    log:
        "logs/gqc/intersect_{category}.log",
    conda:
        "../env/gqc.yml"
    shell:
        "bedtools intersect -a {input.errors} -b {input.annotation} -u > {output} 2>{log}"


# ── Covered-base computation ──────────────────────────────────────────────────

rule gqc_all_repeat_union:
    """Merge all repeat/difficult annotation BEDs into one union for QV computation."""
    input:
        expand(
            config["gqc_error_classification"]["giab_strat_merged"] + "/{strat}.bed",
            strat=["satellite", "segdup", "tandem_repeat", "homopolymer", "mhc", "kir", "vdj"],
        ),
    output:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/all_repeat_union.bed",
    log:
        "logs/gqc/all_repeat_union.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input} \
            | cut -f1-3 \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > {output} 2>{log}
        """


rule gqc_covered_total:
    """Count total covered bases from the benchcovered BED."""
    input:
        "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.benchcovered.HG002v1.1.merged.bed",
    output:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/total.txt",
    log:
        "logs/gqc/covered_total.log",
    shell:
        """
        mkdir -p $(dirname {output})
        awk '!/^track|^#/ {{s += $3 - $2}} END {{print s+0}}' {input} > {output} 2>{log}
        """


rule gqc_covered_cds:
    """Count covered bases that overlap CDS."""
    input:
        covered = "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.benchcovered.HG002v1.1.merged.bed",
        cds     = "data/ref/hg002v1.1.cds_jhu_liftoff.bed",
    output:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/cds.txt",
    log:
        "logs/gqc/covered_cds.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        bedtools intersect -a {input.covered} -b {input.cds} \
            | sort -k1,1 -k2,2n | bedtools merge -i - \
            | awk '{{s += $3 - $2}} END {{print s+0}}' \
            > {output} 2>{log}
        """


rule gqc_covered_repeat:
    """Count covered bases that overlap the all-repeat union."""
    input:
        covered      = "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.benchcovered.HG002v1.1.merged.bed",
        repeat_union = "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/all_repeat_union.bed",
    output:
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/repeat.txt",
    log:
        "logs/gqc/covered_repeat.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        bedtools intersect -a {input.covered} -b {input.repeat_union} \
            | sort -k1,1 -k2,2n | bedtools merge -i - \
            | awk '{{s += $3 - $2}} END {{print s+0}}' \
            > {output} 2>{log}
        """


# ── Classification, QV and HP breakdown ──────────────────────────────────────

rule gqc_error_classification:
    input:
        errors         = "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.errortype.HG002v1.1.bed",
        mnrwithvar     = "analysis_other/GQC/T2T_HG002/assemblybench/T2T_HG002.mononucswithvariants.HG002v1.1.bed",
        covered_total  = "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/total.txt",
        covered_cds    = "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/cds.txt",
        covered_repeat = "analysis_other/GQC/T2T_HG002/gqc_error_classification/covered/repeat.txt",
        flag_cds           = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/cds.bed",
        flag_tandem_repeat = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/tandem_repeat.bed",
        flag_homopolymer   = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/homopolymer.bed",
        flag_satellite     = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/satellite.bed",
        flag_segdup        = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/segdup.bed",
        flag_mhc           = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/mhc.bed",
        flag_kir           = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/kir.bed",
        flag_vdj           = "analysis_other/GQC/T2T_HG002/gqc_error_classification/intersect/vdj.bed",
    output:
        per_error  = "analysis_other/GQC/T2T_HG002/gqc_error_classification/per_error.tsv",
        summary    = "analysis_other/GQC/T2T_HG002/gqc_error_classification/summary.tsv",
        hp_summary = "analysis_other/GQC/T2T_HG002/gqc_error_classification/hp_length_summary.tsv",
    log:
        "logs/gqc/T2T_HG002_error_classification.log",
    conda:
        "../env/gqc.yml"
    run:
        import numpy as np
        import pandas as pd

        log_fh = open(log[0], "w")

        # ── helpers ───────────────────────────────────────────────────────────

        def load_flag_set(bed_path):
            """Return set of (chr, start, end) tuples for pre-intersected errors."""
            hits = set()
            with open(bed_path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("track"):
                        continue
                    f = line.split("\t")
                    hits.add((f[0], int(f[1]), int(f[2])))
            return hits

        def read_int(path):
            return int(open(path).read().strip())

        def parse_alleles(name):
            parts = name.split("_")
            ref = parts[-2]
            alt = parts[-1].split(",")[0]
            return ref, alt, len(alt) - len(ref)

        def qv(n_errors, n_bases):
            if n_bases is None or n_bases <= 0 or n_errors <= 0:
                return np.nan
            return -10 * np.log10(n_errors / n_bases)

        # ── load errors ───────────────────────────────────────────────────────

        print("Loading errors …", file=log_fh)
        err_cols = [
            "chr", "start", "end", "name", "score", "strand",
            "tstart", "tend", "rgb", "error_class", "variant_type", "contig",
        ]
        errors = pd.read_csv(input.errors, sep="\t", header=None, names=err_cols)
        parsed = errors["name"].apply(parse_alleles)
        errors[["ref", "alt", "indel_len"]] = pd.DataFrame(parsed.tolist(), index=errors.index)
        err_keys = list(zip(errors["chr"], errors["start"], errors["end"]))

        # ── homopolymer run lengths from mononucswithvariants ─────────────────

        print("Parsing mononucleotide runs …", file=log_fh)
        mnr_cols = [
            "mnr_chr", "mnr_start", "mnr_end", "mnr_name",
            "s1", "st1", "ts1", "te1", "rgb1",
            "err_chr", "err_start", "err_end", "err_name",
            "s2", "st2", "ts2", "te2", "rgb2",
            "error_class2", "variant_type2", "contig2",
        ]
        mnrv     = pd.read_csv(input.mnrwithvar, sep="\t", header=None, names=mnr_cols)
        mnrv_hit = mnrv[mnrv["err_chr"] != "."].copy()
        mnrv_hit["mnr_length"] = mnrv_hit["mnr_end"] - mnrv_hit["mnr_start"]
        mnrv_hit["mnr_base"]   = mnrv_hit["mnr_name"].apply(lambda n: n.rsplit("_", 1)[-1])
        mnr_lookup = {
            (r.err_chr, r.err_start, r.err_end): (r.mnr_length, r.mnr_base)
            for _, r in mnrv_hit.iterrows()
        }
        errors["in_homopolymer_mnr"] = [k in mnr_lookup for k in err_keys]
        errors["homopolymer_length"] = [mnr_lookup[k][0] if k in mnr_lookup else None for k in err_keys]
        errors["homopolymer_base"]   = [mnr_lookup[k][1] if k in mnr_lookup else None for k in err_keys]

        # ── load pre-computed annotation flags ────────────────────────────────

        print("Loading flag files …", file=log_fh)
        flag_inputs = {
            "cds":           input.flag_cds,
            "tandem_repeat": input.flag_tandem_repeat,
            "homopolymer":   input.flag_homopolymer,
            "satellite":     input.flag_satellite,
            "segdup":        input.flag_segdup,
            "mhc":           input.flag_mhc,
            "kir":           input.flag_kir,
            "vdj":           input.flag_vdj,
        }
        for cat, bed in flag_inputs.items():
            hits = load_flag_set(bed)
            errors[f"in_{cat}"] = [k in hits for k in err_keys]
            print(f"  {cat}: {len(hits):,} overlapping errors", file=log_fh)

        # ── classify (coding > highly_polymorphic > satellite > segdup > TR > HP) ──

        def classify(row):
            if row["in_cds"]:                                    return "coding"
            if row["in_mhc"] or row["in_kir"] or row["in_vdj"]: return "highly_polymorphic"
            if row["in_satellite"]:                              return "satellite"
            if row["in_segdup"]:                                 return "segmental_dup"
            if row["in_tandem_repeat"]:                          return "tandem_repeat"
            if row["in_homopolymer"] or row["in_homopolymer_mnr"]: return "homopolymer"
            return "non_repetitive"

        errors["category"] = errors.apply(classify, axis=1)
        errors.drop(columns=["tstart", "tend", "rgb", "score"]).to_csv(
            output.per_error, sep="\t", index=False
        )
        print(f"Saved {len(errors):,} annotated errors → {output.per_error}", file=log_fh)

        # ── covered-base counts for QV (all / non_repetitive / coding only) ──

        total_cov   = read_int(input.covered_total)
        coding_cov  = read_int(input.covered_cds)
        repeat_cov  = read_int(input.covered_repeat)
        non_rep_cov = total_cov - coding_cov - repeat_cov
        print(f"  total covered:          {total_cov:,} bp", file=log_fh)
        print(f"  coding covered:         {coding_cov:,} bp", file=log_fh)
        print(f"  repeat covered:         {repeat_cov:,} bp", file=log_fh)
        print(f"  non-repetitive covered: {non_rep_cov:,} bp", file=log_fh)

        # ── QV summary table ──────────────────────────────────────────────────

        print("Computing QV …", file=log_fh)
        cat_labels = [
            ("all",                "All regions"),
            ("non_repetitive",     "Non-repetitive"),
            ("coding",             "Coding (CDS)"),
            ("homopolymer",        "Homopolymers"),
            ("tandem_repeat",      "Tandem repeats"),
            ("satellite",          "Satellite repeats"),
            ("segmental_dup",      "Segmental duplications"),
            ("highly_polymorphic", "Highly polymorphic (MHC/KIR/VDJ)"),
        ]
        covered_by_cat = {
            "all": total_cov, "non_repetitive": non_rep_cov, "coding": coding_cov,
            "homopolymer": None, "tandem_repeat": None, "satellite": None,
            "segmental_dup": None, "highly_polymorphic": None,
        }
        rows = []
        for cat, label in cat_labels:
            sub     = errors if cat == "all" else errors[errors["category"] == cat]
            n_bases = covered_by_cat[cat]
            for vtype, vlab in [("all", "All"), ("INDEL", "Indel"), ("SNV", "SNV")]:
                s      = sub if vtype == "all" else sub[sub["variant_type"] == vtype]
                n_all  = len(s)
                n_cons = int((s["error_class"] == "CONSENSUS").sum())
                rows.append(dict(
                    category           = label,
                    variant_type       = vlab,
                    n_bases_covered    = n_bases if n_bases is not None else "NA",
                    n_errors_total     = n_all,
                    pct_of_total       = round(100 * n_all / len(errors), 2) if len(errors) else 0,
                    n_errors_consensus = n_cons,
                    n_errors_phasing   = int(n_all - n_cons),
                    QV_all             = qv(n_all,  n_bases),
                    QV_consensus       = qv(n_cons, n_bases),
                ))
        pd.DataFrame(rows).to_csv(
            output.summary, sep="\t", index=False, float_format="%.2f"
        )

        # ── homopolymer-length breakdown ──────────────────────────────────────

        print("Computing indel counts by homopolymer length …", file=log_fh)
        hp_bins   = [7, 10, 15, 20, 25, 30, 50, 1000]
        hp_labels = ["7-9", "10-14", "15-19", "20-24", "25-29", "30-49", "50+"]
        indels    = errors[errors["variant_type"] == "INDEL"].copy()
        hp_indels = indels[indels["in_homopolymer_mnr"]].copy()
        hp_indels["hp_len_bin"] = pd.cut(
            hp_indels["homopolymer_length"].astype(float),
            bins=hp_bins, labels=hp_labels, right=False,
        )
        hp_rows = []
        for bin_label in hp_labels:
            for ec in ["CONSENSUS", "PHASING", "all"]:
                mask = hp_indels["hp_len_bin"] == bin_label
                if ec != "all":
                    mask = mask & (hp_indels["error_class"] == ec)
                hp_rows.append(dict(hp_length_bin=bin_label, error_class=ec, n_indels=int(mask.sum())))
        non_hp = indels[~indels["in_homopolymer_mnr"]]
        for ec, sub in [
            ("all",       non_hp),
            ("CONSENSUS", non_hp[non_hp["error_class"] == "CONSENSUS"]),
            ("PHASING",   non_hp[non_hp["error_class"] == "PHASING"]),
        ]:
            hp_rows.append(dict(hp_length_bin="non-HP", error_class=ec, n_indels=len(sub)))
        pd.DataFrame(hp_rows).to_csv(output.hp_summary, sep="\t", index=False)

        print("Done.", file=log_fh)
        log_fh.close()
