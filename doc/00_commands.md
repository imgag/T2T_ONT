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

Collect sex informations:
```bash
find assembly/qc/phased_verkko -path '*/T2T*' -name 'sample_sex.txt' -exec sh -c 'echo -n "$1 "; cat "$1"' _ {} \; | sort
```

Test query to update external sample names in NGSD:

```SQL
SELECT
    sample.name_external AS current_name,
    CASE
        WHEN sample.name_external REGEXP '^GE-MED-T2T[0-9]+\\.[0-9]+$' THEN
            REPLACE(SUBSTRING(sample.name_external, 8), '.', '_')
        ELSE
            SUBSTRING(sample.name_external, 8)
    END AS new_name
FROM
    sample
    JOIN processed_sample ON sample.id = processed_sample.sample_id
    JOIN project ON processed_sample.project_id = project.id
WHERE
    project.name = '25006_1422_BEGIN_T2T_GoE'
    AND sample.name_external LIKE 'GE-MED-%';
```

Real queryt to update external sample names in NGSD:

```SQL
UPDATE sample
JOIN processed_sample ON sample.id = processed_sample.sample_id
JOIN project ON processed_sample.project_id = project.id
SET sample.name_external =
    CASE
        -- First replace dots with underscores for patterns like T2T00.1
        WHEN sample.name_external REGEXP '^GE-MED-T2T[0-9]+\\.[0-9]+$' THEN
            REPLACE(
                SUBSTRING(sample.name_external, 8), -- Remove GE-MED- (7 chars)
                '.', -- Replace dots
                '_'  -- With underscores
            )
        -- For all other samples, just remove the GE-MED- prefix
        ELSE
            SUBSTRING(sample.name_external, 8) -- Remove GE-MED- (7 chars)
    END
WHERE
    project.name = '25006_1422_BEGIN_T2T_GoE'
    AND sample.name_external LIKE 'GE-MED-%';
```
