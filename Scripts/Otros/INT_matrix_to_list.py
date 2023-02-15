#!/usr/bin/env python3
import pandas as pd
datos = pd.read_csv("/media/run-projects/Adolfo/MirScience/TCGA/Lung/Immune_INT_matrix_database_genes.tab", sep="\t")


De_data = datos[['gene_id', 'gene_name',"log2FoldChange","padj","DEG"]]
datos = datos.set_index(['gene_id']).drop(columns={'gene_name',"log2FoldChange","padj","DEG"}).drop_duplicates()
De_data.loc[De_data["DEG"]=="Up"]


interactions = pd.DataFrame()
for miRNA in datos.columns:
    for Genes in datos.index:
        if datos.loc[Genes,miRNA] != 0:
            df_in_loop = pd.DataFrame(data={"miRNA":miRNA,"Gene":Genes}, index=[0])
            interactions = interactions.append(df_in_loop)

interactions = interactions.loc[interactions["Gene"].isin(De_data.loc[De_data["DEG"]=="Up"]["gene_id"])]

counts = interactions['miRNA'].value_counts()

interactions = interactions[interactions['miRNA'].isin(counts[counts > 3].index)]
interactions.to_csv("database_miRNAs_UP_genes_list.tab", sep="\t",header=True, index=False)