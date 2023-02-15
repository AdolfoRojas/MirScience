#!/usr/bin/env python3
import pandas as pd 
import numpy as np
IDs = pd.read_csv("Immune_INT_matrix_paper_genes.tab",sep="\t")[["gene_id","DEG","gene_name"]].drop_duplicates() # Immune_genes_Databases_IDs.csv Immune_genes_papers_IDs.csv

INT_matrix = pd.read_csv("../matrix_anotada_full.txt",sep="\t")

Immune_INT_matrix = pd.merge(IDs,INT_matrix,left_on="gene_id",right_on="Gene.stable.ID").drop(columns={"Gene.stable.ID","NCBI.gene..formerly.Entrezgene..ID","HGNC.symbol","gene_name","DEG"}).set_index(['gene_id']).drop_duplicates()
Immune_INT_matrix.columns = Immune_INT_matrix.columns.str.replace(".","-")
sums = Immune_INT_matrix.select_dtypes(pd.np.number).sum().rename('total').to_frame()
Immune_INT_matrix = Immune_INT_matrix[Immune_INT_matrix.columns[Immune_INT_matrix.columns.isin(sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).index)]]
#Immune_INT_matrix.to_csv("Immune_INT_matrix_paper_genes.tab", sep ="\t")
#sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).to_csv("GO_BF_miRNA_INT_counts.tab",sep="\t")

interactions = pd.DataFrame()
for miRNA in Immune_INT_matrix.columns:
    for Genes in Immune_INT_matrix.index:
        if Immune_INT_matrix.loc[Genes,miRNA] != 0:
            df_in_loop = pd.DataFrame(data={"miRNA":miRNA,"Gene":Genes}, index=[0])
            interactions = interactions.append(df_in_loop)


interactions = interactions.merge(IDs, left_on="Gene",right_on="gene_id").drop_duplicates()

IDs2 = pd.read_csv("Immune_INT_matrix_database_genes.tab",sep="\t")[["gene_id","DEG","gene_name"]].drop_duplicates() 

Immune_INT_matrix2 = pd.merge(IDs2,INT_matrix,left_on="gene_id",right_on="Gene.stable.ID").drop(columns={"Gene.stable.ID","NCBI.gene..formerly.Entrezgene..ID","HGNC.symbol","gene_name","DEG"}).set_index(['gene_id']).drop_duplicates()
Immune_INT_matrix2.columns = Immune_INT_matrix2.columns.str.replace(".","-")
sums = Immune_INT_matrix2.select_dtypes(pd.np.number).sum().rename('total').to_frame()
Immune_INT_matrix2 = Immune_INT_matrix2[Immune_INT_matrix2.columns[Immune_INT_matrix2.columns.isin(sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).index)]]
#Immune_INT_matrix.to_csv("Immune_INT_matrix_paper_genes.tab", sep ="\t")
#sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).to_csv("GO_BF_miRNA_INT_counts.tab",sep="\t")

interactions2 = pd.DataFrame()
for miRNA in Immune_INT_matrix2.columns:
    for Genes in Immune_INT_matrix2.index:
        if Immune_INT_matrix2.loc[Genes,miRNA] != 0:
            df_in_loop = pd.DataFrame(data={"miRNA":miRNA,"Gene":Genes}, index=[0])
            interactions2 = interactions2.append(df_in_loop)


interactions2 = interactions2.merge(IDs2, left_on="Gene",right_on="gene_id").drop_duplicates()


###############################################################################################################################################
#                      Generate up and down groups
#Papers
UP_genes_regulators1 = interactions.loc[interactions["DEG"]=="Up"][["miRNA","DEG","gene_name"]].drop_duplicates()
DOWN_genes_regulators1 = interactions.loc[interactions["DEG"]=="Down"][["miRNA","DEG","gene_name"]].drop_duplicates()
#Databases
UP_genes_regulators2 = interactions2.loc[interactions2["DEG"]=="Up"][["miRNA","DEG","gene_name"]].drop_duplicates()
DOWN_genes_regulators2 = interactions2.loc[interactions2["DEG"]=="Down"][["miRNA","DEG","gene_name"]].drop_duplicates()
###############################################################################################################################################
#                      eliminated the genes in common between the UP and down groups 
#Papers
UP_genes_regulators1 = UP_genes_regulators1.loc[~UP_genes_regulators1["gene_name"].isin(UP_genes_regulators2["gene_name"])] 
DOWN_genes_regulators1 = DOWN_genes_regulators1.loc[~DOWN_genes_regulators1["gene_name"].isin(DOWN_genes_regulators2["gene_name"])]
#Databases
UP_genes_regulators2 = UP_genes_regulators2.loc[~UP_genes_regulators2["gene_name"].isin(UP_genes_regulators1["gene_name"])]
DOWN_genes_regulators2 = DOWN_genes_regulators2.loc[~DOWN_genes_regulators2["gene_name"].isin(DOWN_genes_regulators1["gene_name"])]
###############################################################################################################################################
#                      eliminated the microRNAs in common between the UP and Down groups 
#Papers
UP_genes_regulators3 = UP_genes_regulators1.loc[~UP_genes_regulators1["miRNA"].isin(DOWN_genes_regulators1["miRNA"])]
DOWN_genes_regulators3 = DOWN_genes_regulators1.loc[~DOWN_genes_regulators1["miRNA"].isin(UP_genes_regulators1["miRNA"])]
#Databases
UP_genes_regulators4 = UP_genes_regulators2.loc[~UP_genes_regulators2["miRNA"].isin(DOWN_genes_regulators2["miRNA"])]
DOWN_genes_regulators4 = DOWN_genes_regulators2.loc[~DOWN_genes_regulators2["miRNA"].isin(UP_genes_regulators2["miRNA"])]
###############################################################################################################################################
Non_shared_Papers = UP_genes_regulators3.append(DOWN_genes_regulators3)
Non_shared_DB = UP_genes_regulators4.append(DOWN_genes_regulators4)

Data = pd.DataFrame(np.array([["hsa-miR-129-2-3p", len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-129-2-3p"]), len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-129-2-3p"])],["hsa-miR-155-5p",len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-155-5p"]),len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-155-5p"])],["hsa-miR-374a-5p",len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-374a-5p"]),len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-374a-5p"])],["hsa-miR-16-5p",len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-16-5p"]),len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-16-5p"])],["hsa-miR-124-3p",len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-124-3p"]),len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-124-3p"])],["hsa-miR-1-3p",len(Non_shared_Papers.loc[Non_shared_Papers["miRNA"]=="hsa-miR-1-3p"]),len(Non_shared_DB.loc[Non_shared_DB["miRNA"]=="hsa-miR-1-3p"])]]), columns=['miRNA', 'Papers', 'Databases'])





print( + "\t" + str() + " " + str())

###############################################################################################################################################
papers_and_DB = Non_shared_Papers.append(Non_shared_DB)
###############################################################################################################################################
UP_genes_regulators = papers_and_DB.loc[papers_and_DB["DEG"]=="Up"][["miRNA","DEG"]].drop_duplicates()
DOWN_genes_regulators = papers_and_DB.loc[papers_and_DB["DEG"]=="Down"][["miRNA","DEG"]].drop_duplicates()

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Use the venn2 function
venn2(subsets = (len(DOWN_genes_regulators.loc[~DOWN_genes_regulators["miRNA"].isin(UP_genes_regulators["miRNA"])]), 
len(UP_genes_regulators.loc[~UP_genes_regulators["miRNA"].isin(DOWN_genes_regulators["miRNA"])]), 
len(DOWN_genes_regulators.loc[DOWN_genes_regulators["miRNA"].isin(UP_genes_regulators["miRNA"])])
), set_labels = ('Down genes regulator',"Up genes regulator"))
plt.savefig("venn_database_papers_miRNAs.png")

DOWN_genes_regulators.loc[~DOWN_genes_regulators["miRNA"].isin(UP_genes_regulators["miRNA"])]
UP_genes_regulators.loc[~UP_genes_regulators["miRNA"].isin(DOWN_genes_regulators["miRNA"])]

#interactions = interactions.loc[interactions["Gene"].isin(IDs.loc[IDs["DEG"]=="Up"]["gene_id"])]

#counts = interactions['miRNA'].value_counts()

#interactions = interactions[interactions['miRNA'].isin(counts[counts > 3].index)]


