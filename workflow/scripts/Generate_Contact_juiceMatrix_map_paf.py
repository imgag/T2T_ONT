#!/usr/bin/env python
# coding: utf-8

# v1 generate juicebox Medium format contact matrix from map.paf
import sys
import pandas as pd
import numpy as np
import itertools
##from joblib import Parallel, delayed


RvdFfile = sys.argv[1]
Exportfile = sys.argv[2]
print("Generate Contact Matrix for %s"%(RvdFfile) )
colnames = ['read_name', 'read_length', 'read_start', 'read_end', 'strand', 'chrom','chrom_length', 'start', 'end', 'MapQual']
usecols = [0,1,2,3,4,5,6,7,8,11]
RvdF_DF = pd.read_csv(RvdFfile, header=None, index_col=None, sep="\t", usecols=usecols, names=colnames)
RvdF_DF["Position"] = (RvdF_DF["start"].to_numpy() + RvdF_DF["end"].to_numpy())/2
RvdF_DF["Position"] =  RvdF_DF["Position"].astype("int")
RvdF_DF = RvdF_DF.loc[:, ["read_name","chrom", "Position", "strand", "MapQual"]]
selectchrs = [ "chr%d"%i for i in range(1,22+1) ]
selectchrs.extend(["chrX", "chrY"])
RvdF_DF = RvdF_DF[RvdF_DF["chrom"].isin(selectchrs)]
# Change strand
## forward :0;  reverse : -1
RvdF_DF.loc[RvdF_DF.strand=="-", "strand"] = -1
RvdF_DF.loc[RvdF_DF.strand=="+", "strand"] = 0
# ReadGroup Contacts
RvdF_group = RvdF_DF.groupby(by="read_name", as_index=False) 
# Max bins of reads
##Maxbins = RvdF_DF.groupby(by="read_name")["Position"].count().max()
Maxbins = 100
BinDict = {}
for n in range(1, Maxbins):
    lis = [ i for i in range(0, n+1) ]
    BinDict[n+1] = list(itertools.combinations(list(lis),2))

def ReAlignIndex(gDF): 
    gDF.sort_values(by=["chrom", "Position"], ascending=True, inplace=True, ignore_index=True) 
    return (gDF)

def ContactMatrix(gDF): 
    contactDF = pd.DataFrame({"read_name":[], 
                              "str1":[], "chr1":[], "pos1":[],"frag1":[],
                              "str2":[], "chr2":[], "pos2":[],"frag2":[], 
                               "mapq1":[],"mapq2":[]}) 
    if len(gDF) > 1: 
        binidx = BinDict[len(gDF)] 
        i1 = [ t[0] for t in binidx ] 
        i2 = [ t[1] for t in binidx ] 
        contactDF["read_name"] = gDF.iloc[i1]["read_name"].to_list() 
        contactDF["str1"] = gDF.iloc[i1]["strand"].to_numpy("int")
        contactDF["chr1"] = gDF.iloc[i1]["chrom"].to_list() 
        contactDF["pos1"] = gDF.iloc[i1]["Position"].to_numpy("int")
        contactDF["frag1"] = 0
        contactDF["str2"] = gDF.iloc[i2]["strand"].to_numpy("int")
        contactDF["chr2"] = gDF.iloc[i2]["chrom"].to_list() 
        contactDF["pos2"] = gDF.iloc[i2]["Position"].to_numpy("int")
        contactDF["frag2"] = 1
        contactDF["mapq1"] = gDF.iloc[i1]["MapQual"].to_numpy("int")
        contactDF["mapq2"] = gDF.iloc[i2]["MapQual"].to_numpy("int")
    return(contactDF)

def ExportFun(Exportfile, gcontactDF, Nloop):
    # export
    #ContactDF = pd.concat(gcontactDF, axis=0, ignore_index=True)
    ContactDF = pd.DataFrame()
    ContactDF = ContactDF.append(gcontactDF, ignore_index=True) 
    ContactDF = ContactDF.astype({"str1":"int", "str2":"int","pos1":"int","pos2":"int",
                              "frag1":"int","frag2":"int","mapq1":"int","mapq2":"int"})
    ContactDF.to_csv(Exportfile, sep=" ", header=False, index=False, mode="a")
    print("Generated %d pairs contacts in %d reads." %(len(ContactDF), Nloop) )
    return 1

# Generate Matrix and export
gcontactDF = []
Nloop, N, readsNum = 20000, 0, 0
for key, gDF in RvdF_group:
    gcontactDF.append( ContactMatrix( ReAlignIndex(gDF) ) )
    N += 1
    readsNum += 1
    if N % Nloop == 0:
        status = ExportFun(Exportfile, gcontactDF, Nloop)
        gcontactDF = []
        N = 0
        
# The last export
if len(gcontactDF) >= 1:
    status = ExportFun(Exportfile, gcontactDF, N)
    print("%d reads, finished!"%readsNum)


