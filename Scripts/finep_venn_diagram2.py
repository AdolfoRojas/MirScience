python3
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

IDs = pd.read_csv("paper_and_DE_2_databases_genes.tab",sep="\t")[["gene_id","DEG","gene_name"]].drop_duplicates() # Immune_genes_Databases_IDs.csv Immune_genes_papers_IDs.csv

INT_matrix = pd.read_csv("../matrix_anotada_full.txt",sep="\t")

Immune_INT_matrix = pd.merge(IDs,INT_matrix,left_on="gene_id",right_on="Gene.stable.ID").drop(columns={"Gene.stable.ID","NCBI.gene..formerly.Entrezgene..ID","HGNC.symbol","gene_name","DEG"}).set_index(['gene_id']).drop_duplicates()
Immune_INT_matrix.columns = Immune_INT_matrix.columns.str.replace(".","-")
sums = Immune_INT_matrix.select_dtypes(pd.np.number).sum().rename('total').to_frame()
Immune_INT_matrix = Immune_INT_matrix[Immune_INT_matrix.columns[Immune_INT_matrix.columns.isin(sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).index)]]

interactions = pd.DataFrame()
for miRNA in Immune_INT_matrix.columns:
    for Genes in Immune_INT_matrix.index:
        if Immune_INT_matrix.loc[Genes,miRNA] != 0:
            df_in_loop = pd.DataFrame(data={"miRNA":miRNA,"Gene":Genes}, index=[0])
            interactions = interactions.append(df_in_loop)


interactions = interactions.merge(IDs, left_on="Gene",right_on="gene_id").drop(columns={"Gene"}).drop_duplicates()
###############################################################################################################################################
#                      Generate up and down groups
UP_genes_regulators1 = interactions.loc[interactions["DEG"]=="Up"][["miRNA","DEG"]].drop_duplicates()
DOWN_genes_regulators1 = interactions.loc[interactions["DEG"]=="Down"][["miRNA","DEG"]].drop_duplicates()
###############################################################################################################################################
#                      eliminated the microRNAs in common between the UP and Down groups 
UP_genes_regulators2 = UP_genes_regulators1.loc[~UP_genes_regulators1["miRNA"].isin(DOWN_genes_regulators1["miRNA"])]
DOWN_genes_regulators2 = DOWN_genes_regulators1.loc[~DOWN_genes_regulators1["miRNA"].isin(UP_genes_regulators1["miRNA"])]
###############################################################################################################################################
Non_shared = UP_genes_regulators2.append(DOWN_genes_regulators2)
# Use the venn2 function
venn2(subsets = (len(DOWN_genes_regulators1.loc[~DOWN_genes_regulators1["miRNA"].isin(UP_genes_regulators1["miRNA"])]), 
len(UP_genes_regulators1.loc[~UP_genes_regulators1["miRNA"].isin(DOWN_genes_regulators1["miRNA"])]), 
len(DOWN_genes_regulators1.loc[DOWN_genes_regulators1["miRNA"].isin(UP_genes_regulators1["miRNA"])])
), set_labels = ('Down genes regulator',"Up genes regulator"))
plt.savefig("venn_database_papers_miRNAs2.png")