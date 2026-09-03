"""
Assign parent-of-origin labels to assembly contigs using PatMat PofO_Scores.tsv
and a PAF contig-to-chromosome mapping.

The assembler phases each chromosome independently, so HP1 can be maternal on
chr1 and paternal on chr5. For chromosomes where PatMat assigned a parent of
origin, contigs are placed into the correct haplotype:
    maternal → hap1 output  (#1 in PanSN)
    paternal → hap2 output  (#2 in PanSN)
For chromosomes without a PofO assignment (e.g. insufficient iDMR evidence),
the original Verkko haplotype label is preserved (hap1 stays #1, hap2 stays #2).

Outputs (all with PanSN-spec headers):
    hap1     — maternal contigs (assigned) + haplotype1 contigs (unassigned)  {sample}#1#
    hap2     — paternal contigs (assigned) + haplotype2 contigs (unassigned)  {sample}#2#
    combined — hap1 + hap2 + assembly-unassigned contigs from both.fasta      {sample}#0#
    table    — TSV: contig, source_haplotype, mapped_chromosome,
                    pofo_origin, output_haplotype
"""

import argparse
import sys
from collections import defaultdict


def parse_scores(scores_file):
    """Return {chrom: {'HP1': 'maternal'|'paternal', 'HP2': ...}} from PofO_Scores.tsv."""
    assignments = {}
    with open(scores_file) as fh:
        header = fh.readline().rstrip().split("\t")
        try:
            idx_chrom = header.index("Chromosome")
            idx_hp1   = header.index("Origin_HP1")
            idx_hp2   = header.index("Origin_HP2")
        except ValueError as e:
            sys.exit(f"ERROR: unexpected header in {scores_file}: {e}")
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) <= max(idx_chrom, idx_hp1, idx_hp2):
                continue
            chrom = parts[idx_chrom]
            assignments[chrom] = {
                "HP1": parts[idx_hp1].lower(),
                "HP2": parts[idx_hp2].lower(),
            }
    if not assignments:
        sys.exit(f"ERROR: no assignment rows found in {scores_file}")
    return assignments


def parse_paf(paf_file):
    """Return {contig: chromosome} based on largest total alignment block length."""
    aln_len = defaultdict(lambda: defaultdict(int))
    with open(paf_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 11:
                continue
            aln_len[f[0]][f[5]] += int(f[10])
    return {
        contig: max(chroms, key=chroms.get)
        for contig, chroms in aln_len.items()
    }


def iter_fasta(path):
    """Yield (header_line, sequence_lines) for each record in a FASTA."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        yield header, seq


def write_record(fh, header, seq):
    fh.write(header + "\n")
    fh.write("\n".join(seq) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores",   required=True, help="PatMat PofO_Scores.tsv")
    parser.add_argument("--paf",      required=True, help="Contig-to-reference PAF (from both.fasta)")
    parser.add_argument("--hap1",     required=True, help="Verkko haplotype1 FASTA")
    parser.add_argument("--hap2",     required=True, help="Verkko haplotype2 FASTA")
    parser.add_argument("--both",     required=True, help="Verkko assembly.both.fasta (all contigs)")
    parser.add_argument("--out-hap1", required=True, help="Output hap1 FASTA (PanSN #1)")
    parser.add_argument("--out-hap2", required=True, help="Output hap2 FASTA (PanSN #2)")
    parser.add_argument("--combined", required=True, help="Output combined FASTA (PanSN)")
    parser.add_argument("--table",    required=True, help="Output contig assignment TSV")
    parser.add_argument("--sample",   required=True, help="Sample name for PanSN prefix")
    args = parser.parse_args()

    chrom_pofo   = parse_scores(args.scores)
    contig_chrom = parse_paf(args.paf)

    # Infer sex from chrY presence: males have Y contigs in the assembly.
    is_male = any(chrom == "chrY" for chrom in contig_chrom.values())
    print(f"Inferred sex: {'male' if is_male else 'female'}")
    if is_male:
        print("  chrX contigs will be forced to maternal (hap1)")
        print("  chrY contigs will be forced to paternal (hap2)")

    counts = {"hap1": 0, "hap2": 0, "asm_unassigned": 0}
    pofo_counts = {"maternal": 0, "paternal": 0, "unassigned": 0, "sex_chrom": 0}
    seen_contigs = set()

    with open(args.out_hap1, "w") as h1_fh, \
         open(args.out_hap2, "w") as h2_fh, \
         open(args.combined, "w") as comb_fh, \
         open(args.table,    "w") as tbl_fh:

        tbl_fh.write(
            "contig\tsource_haplotype\tmapped_chromosome\tpofo_origin\toutput_haplotype\n"
        )

        # --- hap1 and hap2: PofO assignment or original label ---
        for src_hap, hap_key, fasta_path in [
            ("haplotype1", "HP1", args.hap1),
            ("haplotype2", "HP2", args.hap2),
        ]:
            for header, seq in iter_fasta(fasta_path):
                contig = header[1:].split()[0]
                seen_contigs.add(contig)
                chrom  = contig_chrom.get(contig)

                if is_male and chrom == "chrX":
                    # Single X in males is always maternal
                    pofo_origin = "maternal"
                    out_hap = "hap1"
                    pofo_counts["sex_chrom"] += 1
                elif is_male and chrom == "chrY":
                    # Y is always paternal
                    pofo_origin = "paternal"
                    out_hap = "hap2"
                    pofo_counts["sex_chrom"] += 1
                elif chrom is not None and chrom in chrom_pofo:
                    pofo_origin = chrom_pofo[chrom][hap_key]
                    out_hap = "hap1" if pofo_origin == "maternal" else "hap2"
                else:
                    pofo_origin = "unassigned"
                    out_hap = "hap1" if src_hap == "haplotype1" else "hap2"

                pofo_counts[pofo_origin] += 1
                counts[out_hap] += 1

                pansn_id     = "1" if out_hap == "hap1" else "2"
                pansn_header = f">{args.sample}#{pansn_id}#{contig}"

                dest = h1_fh if out_hap == "hap1" else h2_fh
                write_record(dest,    pansn_header, seq)
                write_record(comb_fh, pansn_header, seq)

                tbl_fh.write(
                    f"{contig}\t{src_hap}\t{chrom or 'NA'}\t"
                    f"{pofo_origin}\t{out_hap}\n"
                )

        # --- assembly-unassigned contigs from both.fasta (not in hap1/hap2) ---
        for header, seq in iter_fasta(args.both):
            contig = header[1:].split()[0]
            if contig in seen_contigs:
                continue
            chrom = contig_chrom.get(contig)
            counts["asm_unassigned"] += 1
            pansn_header = f">{args.sample}#0#{contig}"
            write_record(comb_fh, pansn_header, seq)
            tbl_fh.write(
                f"{contig}\tasm_unassigned\t{chrom or 'NA'}\t"
                f"asm_unassigned\tasm_unassigned\n"
            )

    # --- Stats ---
    total = sum(counts.values())
    print(f"Contig assignment summary ({total} total):")
    print(f"  hap1 (#1) output      : {counts['hap1']}")
    print(f"  hap2 (#2) output      : {counts['hap2']}")
    print(f"    of which maternal   : {pofo_counts['maternal']}")
    print(f"    of which paternal   : {pofo_counts['paternal']}")
    print(f"    of which sex chrom  : {pofo_counts['sex_chrom']} (chrX/Y hardcoded)")
    print(f"    of which unassigned : {pofo_counts['unassigned']} (original haplotype kept)")
    print(f"  asm-unassigned (#0)   : {counts['asm_unassigned']} (combined only)")

    print("\nPer-chromosome PofO assignments:")
    print(f"  {'Chromosome':<16}  {'HP1→':<14}  {'HP2→'}")
    for chrom, vals in sorted(chrom_pofo.items()):
        hp1_out = "hap1 (#1)" if vals["HP1"] == "maternal" else "hap2 (#2)"
        hp2_out = "hap1 (#1)" if vals["HP2"] == "maternal" else "hap2 (#2)"
        print(f"  {chrom:<16}  {vals['HP1']} → {hp1_out:<14}  {vals['HP2']} → {hp2_out}")

    unassigned_chroms = sorted(set(contig_chrom.values()) - set(chrom_pofo.keys()))
    if unassigned_chroms:
        print(f"\nChromosomes without PofO assignment (original labels kept): "
              f"{', '.join(unassigned_chroms)}")


if __name__ == "__main__":
    main()
