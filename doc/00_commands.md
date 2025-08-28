# Bash commands

Collection of instructions to run different parts and tools of the analysis pipelines.

## MinKnow Read QC

To copy run reports:

```bash
find /mnt/storage3/raw_data/MINERVA/25006 -name "report_*.html" -exec cp {} /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/doc/run_reports/ \;
```

Merge reports:

```bash
python workflow/scripts/09_parse_minknow_reports.py doc/run_reports/
```

Export external sample names from NGSD

```bash
python workflow/scripts/19_query_ngsd_for_samplename.py
```

Create plots:

```bash
Rscript workflow/scripts/10_plot_minknow_reports.R doc/run_reports/run_summary.csv doc/tables/flowcell_biological_sample.tsv doc/img/
```