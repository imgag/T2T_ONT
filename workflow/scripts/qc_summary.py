import os, sys
import pandas as pd

path_to_qc_folder = sys.argv[1] #"assembly/qc/phased_verkko"
#usage: python qc_summary.py /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/qc/phased_verkko

os.chdir(path_to_qc_folder)
wp = os.getcwd()
qc_fol = os.listdir(wp)

def get_var(f, lineage):
    run_name = f.replace("published_", "")
    run_name = run_name.replace("_", " ")
    all_stat = {}
    if lineage == "pat":
        lin = "haplotype1"
    else:
        lin = "haplotype2"
    ng50 = None
    nga50 = None
    scaffold_count = None
    longest_scaffold = None
    shortest_scaffold = None

    #getting NG stat from amstat file
    with open(os.path.join(f, f"qc_paftools_asmstat.{lin}.txt"), "r") as file:
        for line in file:
            if line.startswith("NG50"):
                ng50 = line.split("\t")[1]
            if line.startswith("NGA50"):
                nga50 = line.split("\t")[1]
    
    #getting scaffold stats from scaffolds length stats file
    with open(os.path.join(f, f"scaffold_lengths.{lin}.txt"), "r") as file:
        length=[]
        for line in file:
            if line.startswith("haplotype"):
                #print(line)
                length.append(int(line.split("\t")[1]))
        scaffold_count = len(length)
        longest_scaffold = max(length)
        shortest_scaffold = min(length)
        total_length = sum(length)

    #getting reference comparision with asmgene output
    with open(os.path.join(f, f"qc_paftools_asmgene.{lin}.txt"), "r") as file:
        sc_gene_lost=None
        dup_gene=None
        mc_gene_lost=None
        for line in file:
            if line.startswith("X\tfull_sgl"):
                total_single_copy=int(line.split("\t")[3])
                sc_gene_lost=(( (int(line.split("\t")[2]) - int(line.split("\t")[3])) / int(line.split("\t")[2])))
                sc_gene_lost = round(sc_gene_lost, 3)
            if line.startswith("X\tfull_dup"):
                dup_gene=int(line.split("\t")[3]) / total_single_copy
                dup_gene=round(dup_gene,3)
            if line.startswith("X\tdup_cnt"):
                mc_gene_lost= 1 - (int(line.split("\t")[3])) / int(line.split("\t")[2])
                mc_gene_lost=round(mc_gene_lost,3)
        asmgene=[sc_gene_lost, dup_gene, mc_gene_lost]

    #NEED TO RUN MERQURY FIRST

    # with open(os.path.join(f, f"QV_score.qv"), "r") as file:
    #     for line in file:
    #         if line.startswith(f"assembly.{lin}"):
    #             qv = line.split("\t")[3]
    #             qv= round(float(qv), 3)

    #getting whatshap comparision stats
    with open(os.path.join(f, "whatshap_compare.tsv"), "r") as file:
        line_count = sum(1 for line in file)
    with open(os.path.join(f, "whatshap_compare.tsv"), "r") as file:
        hamming_error_rate = None
        switch_error_rate = None
        #line_count = sum(1 for line in file)
        if line_count == 2:
            for line in file:
                if line.startswith("SAMPLE"):
                    hamming_error_rate = float(line.split("\t")[14])
                    switch_error_rate = float(line.split("\t")[10])
                    print(hamming_error_rate)
                    print(switch_error_rate)
        if line_count > 2:
            hamms = []
            switchs = []
            for line in file:
                if line.startswith("SAMPLE"):
                    hamms.append(float(line.split("\t")[14]))
                    switchs.append(float(line.split("\t")[10]))    
            hamming_error_rate = sum(hamms) / len(hamms)
            switch_error_rate = sum(switchs) / len(switchs)    

    all_stat["Run"] = run_name
    all_stat["Phase"] = lineage            
    all_stat["NG50"] = int(ng50)
    all_stat["NGA50"] = int(nga50)
    all_stat["scaffold_count"] = scaffold_count
    all_stat["longest_scaffold"] = longest_scaffold
    all_stat["shortest_scaffold"] = shortest_scaffold
    all_stat["total_length"] = total_length

    #need to run merqury first
    #all_stat["QV"] = qv
    
    all_stat["asmgene"] = asmgene
    return all_stat



dict_list = []
for folds in qc_fol:
    dict_list.append(get_var(folds, "pat"))
    dict_list.append(get_var(folds, "mat"))

df=pd.DataFrame(dict_list)
df.to_csv(os.path.join(path_to_qc_folder,"qc_summary_all.tsv"), sep="\t", index=False)