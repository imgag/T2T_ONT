def get_ref_genome(wc):
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
                    "done": f"assembly/output/verkko/{wc['asm']}/create_porec_scaffold.done",
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