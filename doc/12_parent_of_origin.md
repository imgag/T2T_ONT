# Parent-of-Origin Assignment

Determines maternal/paternal origin of assembly contigs using CpG methylation
at imprinted DMRs via [PatMat](https://github.com/vahidAK/PatMat).

## Rules

```
patmat_filter_vcf     → filter VCF to diploid chromosomes (autosomes ± chrX by sex)
patmat                → run PatMat, bgzip+tabix VCF output
assign_haplotype_pofo → per-contig origin assignment, write final FASTAs + table
```

## Outputs

| File | Description |
|------|-------------|
| `{sample}_PofO_Scores.tsv` | Per-chromosome assignment scores |
| `{sample}_PofO_Assignment.vcf.gz` | Variants with PofO tags |
| `{sample}_PofO_Tagged.cram` | Reads tagged with maternal/paternal origin |
| `assembly.hap1.pansn.fasta` | Maternal (assigned) / haplotype1 (unassigned) — PanSN `#1` |
| `assembly.hap2.pansn.fasta` | Paternal (assigned) / haplotype2 (unassigned) — PanSN `#2` |
| `assembly.combined.pansn.fasta` | All contigs: hap1 + hap2 + assembly-unassigned (`#0`) |
| `contig_pofo_assignment.tsv` | Per-contig: source haplotype, chromosome, PofO origin, output haplotype |

PanSN headers: `{sample}#1#contig`, `{sample}#2#contig`, `{sample}#0#contig`.

## Known limitations

**chr17** cannot be assigned. The chromosome has only one iDMR in the
reference panel (CHRNE/C17orf107, chr17:4,790,606–4,792,558) and it fails the
minimum CpG difference threshold (observed differential methylation ~19% vs
required 50%). Chr17 contigs are retained with their original Verkko haplotype
labels in the final assembly.
