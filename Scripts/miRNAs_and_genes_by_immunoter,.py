#!/usr/bin/env python3
import pandas as pd 
gmt = pd.read_csv("immuno_GO_BP.gmt", sep = ")", header = None)
gmt[0] = gmt[0]+")"
gmt[1] = gmt[1].str.replace("\t\t","").str.replace("\t", "; ").str.strip("; ")
IDs = pd.read_csv("TCGA/Lung/2_gene_ID_to_gene_symbol_TCGA.tab",sep="\t")
interactions = pd.read_csv("TCGA/Lung/GO_BP_miRNAs_UP_genes_list.tab",sep="\t")
interactions = pd.merge(interactions,IDs, left_on="Gene",right_on="gene_id").drop_duplicates()

df_attr = dict()
df_attr["Term"] = []
df_attr["Genes"] = []
df_attr["miRNAs"] = []
df_attr["N°Genes"] = []
df_attr["N°miRNAs"] = []
for Index in range(0,len(gmt[0])):
    print(Index)
    term = gmt[0][Index]
    genes_UP_in_term = []
    for gene in gmt[1][Index].split("; "):
        if gene in interactions["gene_name"].values:
            genes_UP_in_term.append(gene)
            print(genes_UP_in_term)
    miRNAs_with_interaction = interactions.loc[interactions["gene_name"].isin(genes_UP_in_term)]["miRNA"].unique().tolist()
    n_genes = len(genes_UP_in_term)
    n_miRNAs = len(miRNAs_with_interaction)
    df_attr["Term"].append(term)
    df_attr["Genes"].append(genes_UP_in_term)
    df_attr["miRNAs"].append(miRNAs_with_interaction)
    df_attr["N°Genes"].append(n_genes)
    df_attr["N°miRNAs"].append(n_miRNAs)

df_attr_df = pd.DataFrame(df_attr)
df_attr_df.to_csv("summary_immuno_GO_BP_terms_miRNA_interactions.tsv",sep = "\t", index = False)

miRNAs_attr = dict()
miRNAs_attr["miRNA"] = []
miRNAs_attr["Genes"] = []
miRNAs_attr["N°Genes"] = []
miRNAs_attr["Immuno_pathways_selected"] = []
miRNAs_attr["N°Immuno_pathways_selected"] = []
for miRNA in interactions["miRNA"].unique():
    genes_with_interaction = interactions.loc[interactions["miRNA"] == miRNA]["gene_name"].unique().tolist()
    n_genes_with_interaction = len(genes_with_interaction)
    immunopaths = df_attr_df.loc[df_attr_df["miRNAs"].astype(str).str.contains(miRNA)]["Term"].unique().tolist()
    n_immunopaths = len(immunopaths)
    miRNAs_attr["miRNA"].append(miRNA)
    miRNAs_attr["Genes"].append(genes_with_interaction)
    miRNAs_attr["N°Genes"].append(n_genes_with_interaction)
    miRNAs_attr["Immuno_pathways_selected"].append(immunopaths)
    miRNAs_attr["N°Immuno_pathways_selected"].append(n_immunopaths)

df_miRNAs_attr = pd.DataFrame(miRNAs_attr)
df_miRNAs_attr.to_csv("summary_miRNA_interactions.tsv",sep = "\t", index = False)

interactions = interactions.drop(columns="Gene")
interactions_selected = interactions.loc[(interactions["miRNA"] == "hsa-miR-129-2-3p")|(interactions["miRNA"] == "hsa-miR-374a-5p")|(interactions["miRNA"] == "hsa-miR-124-3p")|(interactions["miRNA"] == "hsa-miR-155-5p")|(interactions["miRNA"] == "hsa-miR-16-5p")|(interactions["miRNA"] == "hsa-miR-1-3p")]

interactions_selected.to_csv("selected_miRNA_interactions.tsv",sep = "\t", index = False)