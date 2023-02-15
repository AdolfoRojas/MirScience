#!/usr/bin/env python3
from bioinfokit import analys, visuz
import pandas as pd
# load dataset as pandas dataframe
file = "MirScience/TCGA/Pancreatic/mRNAs_TCGA-PAAD_Normal_vs_Tumoral_DE.tab"
df = pd.read_csv(file)
df["padj"] = df["padj"].fillna(1)
visuz.GeneExpression.volcano(df=df, lfc='log2FoldChange', pv='padj', color=("#E10600FF", "grey", "#00239CFF"),r = 500,figname = file, figtype='png')