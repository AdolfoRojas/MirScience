#!/usr/bin/env python3
#%%
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np
import pandas as pd
# Data info
Disease = "nash"
# Project_MethodType Organism_Name
# taxon gdsType PubMedIds
# is_public Library_Strategy Library_Selection* Library_Layout Organism_name Platform
#datos = pd.read_csv("/media/run-projects/Adolfo/MirScience/Datasets/BioProject/"+Disease+"v2.tsv",sep="\t", header= 0)
#datos = pd.read_csv("/media/run-projects/Adolfo/MirScience/Datasets/BioProject/"+Disease+"_GEO_projects.tsv",sep="\t", header= 0)
datos = pd.read_csv("/media/run-projects/Adolfo/MirScience/Datasets/BioProject/"+Disease+"_SRA_Samples.tsv",sep="\t", header= 0)
datos1 = datos["Platform"].value_counts().to_frame()
datos3 = datos["Platform"].value_counts(normalize=True).to_frame()
fig = plt.figure(figsize=(9, 9))
def absolute_value(val):
    a  = np.round(val/100.*ratios.sum(), 0)
    a = int(a)
    return a

ratios = datos1["Platform"].values#.tolist()
ratios2 = datos3["Platform"].values#.tolist()
labels1 = datos1.index.str.replace("_"," ").values.tolist()
angle = -180 * ratios2[0]+180
cmap = plt.get_cmap("tab20")
plt.pie(ratios, autopct=absolute_value, startangle=angle, # autopct='%1.1f%%'
        labels=None, textprops={'fontsize': 18}, colors =cmap(np.arange(0,len(labels1))))

first_legend  =fig.legend(labels1, loc="lower right")
plt.title("Platform", fontsize=40)
plt.savefig('/media/run-projects/Adolfo/MirScience/Datasets/BioProject/Pie_plot_Projects9.png', dpi = 300)

# %%
