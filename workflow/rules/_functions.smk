def get_ref_genome(wc):
# Legacy function that was used when we assembled only a single chromosome.
# Returns the ref genome in other cases.
    import re
    ref = config["ref"]
    # print(ref)
    match = re.search(r"chr\d+", str(wc.asm))
    # print(match)
    if match:
        # print(match)
        prefix = config["ref"].replace("fasta", "")
        ref = f"{prefix}{match.group(0)}.fasta"
    return ref


def get_assembly_output(wc):
    # Treat undefined isphased as "phased"
    isphased = wc.get("isphased", "phased")

    # Not defined for variant calling

    if isphased == "unphased":
        if wc["tool"] == "verkko":
            return {
                "assembly": f"assembly/output/verkko_unphased/{wc['asm']}/assembly.fasta"
            }
        else:
            raise ValueError(f"Invalid tool for unphased assembly: {wc['tool']}")

    elif isphased == "phased":
        if wc["tool"] == "verkko":
            if not wc["hp"] or wc["hp"] == "unphased" or wc["hp"] == "both":
                return {
                    "assembly": f"assembly/output/verkko/{wc['asm']}/assembly.fasta",
                    "done": f"assembly/output/verkko/{wc['asm']}/create_scaffold.done",
                }
            elif wc["hp"] == "haplotype1":
                return {
                    "assembly": f"assembly/output/verkko/{wc['asm']}/assembly.haplotype1.fasta"
                }
            elif wc["hp"] == "haplotype2":
                return {
                    "assembly": f"assembly/output/verkko/{wc['asm']}/assembly.haplotype2.fasta"
                }
        elif wc["tool"] == "gfase":
            if wc["hp"] == "haplotype1":
                return {
                    "assembly": f"assembly/output/gfase/{wc['asm']}/gfase/phase_0.fasta"
                }
            elif wc["hp"] == "haplotype2":
                return {
                    "assembly": f"assembly/output/gfase/{wc['asm']}/gfase/phase_1.fasta"
                }
        elif wc["tool"] == "hifiasm":
            if not wc["hp"] or wc["hp"] == "unphased" or wc["hp"] == "both":
                return {
                    "assembly": f"assembly/output/hifiasm/{wc['asm']}/assembly.fasta"
                }
            elif wc["hp"] == "haplotype1":
                return {
                    "assembly": f"assembly/output/hifiasm/{wc['asm']}/assembly.haplotype1.fasta"
                }
            elif wc["hp"] == "haplotype2":
                return {
                    "assembly": f"assembly/output/hifiasm/{wc['asm']}/assembly.haplotype2.fasta"
                }
        elif wc["tool"] == "repair":
            if not wc["hp"] or wc["hp"] == "unphased" or wc["hp"] == "both":
                return {
                    "assembly": f"analysis_other/assembly_repairer/{wc['asm']}/repaired_assembly.final.fa"
                }
            elif wc["hp"] == "haplotype1":
                raise ValueError(f"Currently no distinct haplotypes for repaired assembly")
                return {
                    "assembly": []
                }
            elif wc["hp"] == "haplotype2":
                raise ValueError(f"Currently no distinct haplotypes for repaired assembly")
                return {
                    "assembly": []
                }
        elif wc["tool"] == "hprc":
            if wc["hp"] == "haplotype1" or wc["hp"] == "haplotype2":
                return {
                    "assembly": f"data/ref/hprc/{wc['asm']}/{wc['asm']}.{wc['hp']}.fasta"
                }
            else:
                raise ValueError(f"Invalid haplotype {wc['hp']} for hprc assembliues")
        elif wc["tool"] == "apkpolish":
            if  not wc["hp"] or wc["hp"] in ["unphased", "both"]:
                return {
                    "assembly": f"assembly/output/apk_polish/{wc['asm']}/assembly.fasta"
                }
            elif wc["hp"] in ["haplotype1", "haplotype2"]:
                #raise ValueError(f"Currently no distinct haplotypes for polished assembly")
                return {
                    "assembly": []
                }
            else:
                raise ValueError(f"Invalid haplotype {wc['hp']} for apk-polished assembliues")
        else:
            raise ValueError(f"Invalid tool: {wc['tool']}")


def get_assembly_graph_output(wc):
    if wc["tool"] == "verkko":
        return f"assembly/output/verkko/{wc['asm']}/assembly.homopolymer-compressed.noseq.gfa"
    elif wc["tool"] == "gfase":
        return f"assembly/output/gfase/{wc['asm']}/gfase/chained.gfa"


def get_assembly_colors(wc):
    if wc["isphased"] == "phased":
        if wc["tool"] == "verkko":
            return f"assembly/output/verkko/{wc['asm']}/assembly.colors.csv"
        if wc["tool"] == "gfase":
            return f"assembly/output/gfase/{wc['asm']}/gfase/phases.csv"
    elif wc["isphased"] == "unphased":
        return f"assembly/qc/unphased_{wc['tool']}/{wc['asm']}/colors.tsv"
    else:
        raise ValueError(f"Invalid  phasing value: {wc['isphased']}")

def get_chr_list_for_asm(wc):

    file = f"assembly/qc/phased_verkko/{wc.asm}/sample_sex.txt"

    if open(file).read().strip() == "male":
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X", "Y"]]
    else:
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X"]]

    return(chroms)

def get_chrom_size(chrom_sizes_file, chrom):
    """Get chromosome size from chrom.sizes file."""
    with open(chrom_sizes_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] == chrom:
                return int(parts[1])
    raise ValueError(f"Chromosome {chrom} not found in {chrom_sizes_file}")

