#!/usr/bin/env python3
#%%
from bioinfokit import analys, visuz
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DE = pd.read_csv("/media/run-projects/Adolfo/MirScience/TCGA/Lung/DE_all/miRNAs_TCGA-LUAD_Normal_vs_Tumoral_DE_NS.csv",sep=",").rename(columns={"Unnamed: 0":"miRNA"})
DE["log2FoldChange2"] = DE["log2FoldChange"].abs()

datos = DE.copy()
datos = datos.loc[datos["miRNA"].str.contains("hsa-mir-29b")]
datos = datos.set_index(datos["miRNA"]).drop(columns={"miRNA"})


materials = datos.index
x_pos = np.arange(len(materials))
CTEs = datos["log2FoldChange"]
error = datos["lfcSE"]
fig = plt.figure(figsize=(32, 9),dpi=500)
fig, ax = plt.subplots()
ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=3,color = "#c2272d")
ax.set_ylabel('Fold Expression')
ax.set_xticks(x_pos)
ax.set_xticklabels(materials)
plt.rc('xtick', labelsize=4)
ax.set_title('Differentially expressed genes')
plt.xticks(rotation = 90)

# Save the figure and show
plt.tight_layout()
plt.savefig('bar_plot_with_error_bars.png')
plt.show()
#%%
plt.savefig('/media/run-projects/Adolfo/MirScience/TCGA/Lung/DE_all/Bar_plot_DEGs.png',dpi=500)