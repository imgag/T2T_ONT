# Methods

Documentation of analysis methods, tools used in the pipeline, and references for key computational approaches.

## HERRO error correction for assembly input

From the Verkko developers on best practice for HERRO-corrected reads as Verkko input ([source](https://github.com/marbl/verkko/issues/322)):

> You can see some info on the README as well as #321. We usually use longest 50x or whatever can be corrected in a reasonable time with Herro. The suggestions pertaining to minimum k-mer coverage in #321 can be used to address the increasing complexity in the graph with high coverage. There are uncorrected regions with systematic errors that start confirming each other and need to be filtered. You can also try the hifiasm correction as an alternative to herro.
> You should not input corrected reads as UL data. The correction can trim and lose parts of the sequence so you want to allow verkko to fill those gaps with the original UL data. Read length is fine, longer reads will build a more resolved initial graph.

## Assembly polishing

Polishing strategies and resources: https://github.com/arangrhie/T2T-Polish/tree/master

## RepeatMasker

RepeatMasker is configured to use HMMER3.1 with the DFAM database as search engine. The DFAM database (mammal/human) is downloaded from https://www.dfam.org/releases/current/families/FamDB/ in H5 format.

## HG002 cell line culturing

The HG002 cell line was thawed and cultured according to the manufacturer's protocol in RPMI (Sigma Aldrich) supplemented with 2 mM L-Glutamine, 100 U/mL penicillin, 100 µg/mL streptomycin, and 10% heat-inactivated FCS (Thermo Fisher Scientific) at 37 °C with 5% CO2.

## GQC error classification

Assembly errors identified by GenomeQC (GQC) against the HG002v1.1 benchmark were annotated using the GIAB v4.2 genome stratifications for HG002 [1]. Because GQC reports errors in HG002v1.1 haplotype coordinates (chr\*_MATERNAL / chr\*_PATERNAL), the haplotype-specific stratification BED files were used directly rather than a reference-genome liftover. Maternal and paternal BED files for each stratum were merged with `bedtools merge` before intersection. Coding sequence annotations were derived from the JHU Liftoff v0.6 RefSeq liftover onto HG002v1.1. rDNA loci were extracted from the cenSat v2.0 annotation and are reported for completeness, but are subsumed under the satellite category for all summary statistics because they represent a structurally distinct subset of satellite arrays rather than an independent error class.

Each error was assigned to exactly one category according to the following priority order: coding > highly polymorphic (MHC, KIR, VDJ) > satellite > segmental duplication > tandem repeat > homopolymer > non-repetitive. The hierarchy reflects both functional importance and the reliability of the benchmark in each context. Coding regions are placed first because errors there carry direct functional consequence regardless of their repeat context. Highly polymorphic loci (MHC, KIR, VDJ) follow immediately because extreme allelic diversity in these regions means that apparent errors may partly reflect genuine population variation rather than assembly inaccuracy, and their benchmark accuracy is accordingly lower. Satellite and segmental duplication regions are ranked next as structurally the most challenging for both assembly and benchmarking. Tandem repeats precede homopolymers because homopolymers are mechanistically a specific subtype of tandem repeat (single-base runs) with a distinct nanopore sequencing error mode.

Assembly quality values (QV = −10 log₁₀(e/b), where e is the error count and b the number of covered benchmark bases) were computed only for the all-errors, non-repetitive, and coding categories. QV is not reported for repeat categories because the covered-base denominator is unreliable in those contexts: alignment ambiguity causes systematic under-mapping of reads to repeat arrays, GQC explicitly excludes many repeat-rich regions from the benchmark callable set, and the benchmark itself has reduced sensitivity in satellite and highly polymorphic loci. The non-repetitive covered-base count was derived as total benchmark coverage minus the union of CDS and all repeat strata.

[1] Dwarshuis, N., Kalra, D., McDaniel, J. et al. The GIAB genomic stratifications resource for human reference genomes. Nat Commun 15, 9029 (2024). https://doi.org/10.1038/s41467-024-53260-y

## Assembly issue repeat-context classification

Assembly issue regions called by NucFlag (COLLAPSE, COLLAPSE_OTHER, COLLAPSE_VAR, HET, MISJOIN) and Flagger (Col, Dup, Err, Hap) were classified into three mutually exclusive repeat context categories using RepeatMasker annotations in Verkko assembly coordinates (no liftover required):

1. **Satellite** — regions overlapping any `Satellite/*` RepeatMasker entry (alpha-satellite, HSAT, acro-satellite, subtelomeric repeats). BED derived by filtering the RepeatMasker `.out` file on column 11 matching `^Satellite`.
2. **Other repeat** — regions overlapping any RepeatMasker-annotated element after subtracting satellite entries (captures LINE/L1, SINE/Alu, LTR, DNA transposons, simple repeats, etc.).
3. **Non-repetitive** — regions with no RepeatMasker annotation (subtracted all annotated repeats).

Priority was applied in that order (satellite takes precedence over other repeat; other repeat takes precedence over non-repetitive). Regions on unplaced/unassigned contigs were excluded; only primary haplotype contigs (haplotype1-\*, haplotype2-\*) were considered. Regions < 500 bp were excluded from the per-region size analysis but not from the summary counts.
