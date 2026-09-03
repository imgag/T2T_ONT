python workflow/scripts/14_rfhap_preprocess_diploid_assembly.py \
    --hap1 assembly/qc/phased_verkko/TUE_02_03UL/haplotype1.fa \
    --hap2 assembly/qc/phased_verkko/TUE_02_03UL/haplotype2.fa \
    --prefix1 hap1_ \
    --prefix2 hap2_ \
    --output analysis_other/rfhap/04_3UL_bin_asm/TUE_02_03UL.fragmented

./bin/nextflow run bin/rfhap-1.0_mod/rfhap.nf \
    --paternal_reads analysis_other/rfhap/04_3UL_bin_asm/paternal.txt \
    --maternal_reads analysis_other/rfhap/04_3UL_bin_asm/maternal.txt \
    --child_reads analysis_other/rfhap/04_3UL_bin_asm/child.txt \
    --outdir analysis_other/rfhap/04_3UL_bin_asm

python workflow/scripts/15_rfhap_merge_phased_chunks.py \
  -a analysis_other/rfhap/04_3UL_bin_asm/rf/haplotypes/TUE_02_03UL.fragmented_fastkm_matrix.hapA.txt \
  -b analysis_other/rfhap/04_3UL_bin_asm/rf/haplotypes/TUE_02_03UL.fragmented_fastkm_matrix.hapB.txt \
  -u analysis_other/rfhap/04_3UL_bin_asm/rf/haplotypes/TUE_02_03UL.fragmented_fastkm_matrix.hapU.txt \
  -g analysis_other/rfhap/04_3UL_bin_asm/TUE_02_03UL.fragmented.gfa \
  -f analysis_other/rfhap/04_3UL_bin_asm/TUE_02_03UL.fragmented.fasta \
  -o analysis_other/rfhap/04_3UL_bin_asm/TUE_02_03UL.fragmented.merged \
  --bridge_unassigned \
  --max_unassigned_bridge 3