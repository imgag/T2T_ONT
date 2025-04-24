./bin/nextflow run bin/rfhap-1.0_mod/rfhap.nf \
    --paternal_reads analysis_other/rfhap/03_3UL_mod/paternal.txt \
    --maternal_reads analysis_other/rfhap/03_3UL_mod/maternal.txt \
    --child_reads analysis_other/rfhap/03_3UL_mod/child.txt \
    --modbam analysis_other/rfhap/03_3UL_mod/modbam.txt \
    --outdir analysis_other/rfhap/03_3UL_mod