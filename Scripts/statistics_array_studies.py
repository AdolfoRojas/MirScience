#!/usr/bin/env python3
import os
import sys
from os import walk
import pandas as pd

Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
files = Output_dir_files.loc[Output_dir_files.file.str.contains("_dif_exp.tsv")].copy()


Summary = pd.DataFrame()
for file in files.file:
    name = file.replace("_dif_exp.tsv","")
    datos = pd.read_csv(name+"/"+name+"_dif_exp_NS.tsv",sep="\t")[["ID","logFC","AveExpr","t","P.Value","adj.P.Val","B","diffexpressed"]]    
    N_genes = len(datos["ID"])
    N_DEG = len(datos.loc[datos["adj.P.Val"] < 0.05]["ID"])
    N_UP = len(datos.loc[(datos["adj.P.Val"] < 0.05) & (datos["logFC"] > 0)]["ID"])
    N_DOWN = len(datos.loc[(datos["adj.P.Val"] < 0.05) & (datos["logFC"] < 0)]["ID"])
    df_in_loop = pd.DataFrame(data={"Study":name,"Total genes":N_genes,"Number of DEGs":N_DEG, "UP-regulated genes":N_UP, "DOWN-regulated genes":N_DOWN}, index=[0])
    Summary = Summary.append(df_in_loop)

Summary.to_csv("Summary_Array_estudies.tab", sep="\t",header=True, index=False)