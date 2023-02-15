#!/usr/bin/env python3
import os
import sys
from os import walk
import pandas as pd

Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
files = Output_dir_files.loc[Output_dir_files.file.str.contains("_Control_vs_Case_DE_NS.csv")].copy()


Summary = pd.DataFrame()
for file in files.file:
    datos = pd.read_csv(file)
    name = file.replace("_Control_vs_Case_DE_NS.csv","")
    N_genes = len(datos["gene_name"])
    N_DEG = len(datos.loc[datos["padj"] < 0.05]["gene_name"])
    N_UP = len(datos.loc[(datos["padj"] < 0.05) & (datos["log2FoldChange"] > 0)]["gene_name"])
    N_DOWN = len(datos.loc[(datos["padj"] < 0.05) & (datos["log2FoldChange"] < 0)]["gene_name"])
    df_in_loop = pd.DataFrame(data={"Study":name,"Total genes":N_genes,"Number of DEGs":N_DEG, "UP-regulated genes":N_UP, "DOWN-regulated genes":N_DOWN}, index=[0])
    Summary = Summary.append(df_in_loop)

Summary.to_csv("Summary_RNAseq_estudies.tab", sep="\t",header=True, index=False)
