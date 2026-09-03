import pandas as pd

# Read file list — source_path and alias derived per category
_df = pd.read_csv("data/upload_folder/GHGA/ghga_upload_file_list.tsv", sep="\t")
_df = _df[_df["found"] == True]

def _upload_name(row):
    if row["category"] == "Raw BAM":
        return row["filename"]
    elif row["category"] == "Assembly FASTA":
        return f"{row['sample']}.assembly.{row['type']}.fasta"
    elif row["category"] == "Phased VCF":
        return f"{row['sample']}.phased.vcf.gz"

def _alias(row):
    if row["category"] == "Raw BAM":
        run_id = row["filename"].split(".")[0].split("_")[-1]
        return f"file_{run_id}"
    elif row["category"] == "Assembly FASTA":
        return f"pfile_{row['sample']}_{row['type']}"
    elif row["category"] == "Phased VCF":
        return f"pfile_{row['sample']}_variants"

_source_map = {_upload_name(r): r["source_path"] for _, r in _df.iterrows()}
_alias_map  = {_upload_name(r): _alias(r)          for _, r in _df.iterrows()}

UPLOAD_NAMES = list(_source_map.keys())

wildcard_constraints:
    name = "[^/]+"

# Before running, export your GHGA access token:
#   export GHGA_TOKEN="<token from portal>"
# Key paths and optional proxy are set in workflow/config.yml

rule all_ghga:
    input:
        expand("data/upload_folder/GHGA/uploaded/{name}.done", name=UPLOAD_NAMES)


rule ghga_upload:
    input:
        lambda wc: _source_map[wc.name]
    output:
        touch("data/upload_folder/GHGA/uploaded/{name}.done")
    params:
        alias   = lambda wc: _alias_map[wc.name],
        pub_key = config["ghga_pub_key"],
        sec_key = config["ghga_sec_key"],
        proxy   = "HTTPS_PROXY={p} HTTP_PROXY={p} ALL_PROXY={p}".format(p=config["ghga_proxy"]) if config.get("ghga_proxy") else ""
    log:
        "logs/ghga_upload/{name}.log"
    threads:
        1
    shell:
        """
        tmp=$(mktemp --suffix=.tsv)
        echo -e "$(realpath {input})\t{params.alias}" > $tmp
        echo "$GHGA_TOKEN" | {params.proxy} uvx --with socksio ghga-connector batch-upload \
            --tsv $tmp \
            --my-public-key-path {params.pub_key} \
            --my-private-key-path {params.sec_key} \
            >> {log} 2>&1
        rm $tmp
        grep -qE "All files uploaded successfully|All files are already uploaded" {log} || {{ echo "Upload failed for {params.alias}, check {log}"; exit 1; }}
        """
