#!/usr/bin/env python3
import pandas as pd 
IDs = pd.read_csv("Immune_genes_papers_IDs.csv")
#gencode_eq = pd.read_csv("2_gene_ID_to_gene_symbol_TCGA.tab", sep="\t")
#IDs = gencode_eq.loc[gencode_eq["gene_name"].isin(IDs["gene_name"])]
#IDs.to_csv("Immune_genes_Databases_IDs.csv")
INT_matrix = pd.read_csv("../matrix_anotada_full.txt",sep="\t")
DE_info = pd.read_csv("mRNAs_TCGA-LUAD_Normal_vs_Tumoral_DE_NS.csv")[["gene_id","gene_name","log2FoldChange","padj"]]
DE_info = DE_info.loc[DE_info["gene_id"].isin(IDs["gene_id"])]
DE_info["DEG"] = "NS"
DE_info.loc[(DE_info["padj"] < 0.05) & (DE_info["log2FoldChange"] < 0), "DEG"] = "Down"
DE_info.loc[(DE_info["padj"] < 0.05) & (DE_info["log2FoldChange"] > 0), "DEG"] = "Up"
Immune_INT_matrix = pd.merge(DE_info,INT_matrix,left_on="gene_id",right_on="Gene.stable.ID").drop(columns={"Gene.stable.ID","NCBI.gene..formerly.Entrezgene..ID","HGNC.symbol"}).set_index(['gene_id', 'gene_name',"log2FoldChange","padj","DEG"]).drop_duplicates()
Immune_INT_matrix.columns = Immune_INT_matrix.columns.str.replace(".","-")
sums = Immune_INT_matrix.select_dtypes(pd.np.number).sum().rename('total').to_frame()
Immune_INT_matrix = Immune_INT_matrix[Immune_INT_matrix.columns[Immune_INT_matrix.columns.isin(sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).index)]]
Immune_INT_matrix.to_csv("Immune_INT_matrix_paper_genes.tab", sep ="\t")
sums.loc[sums["total"] > 0].sort_values(by="total", ascending= False).to_csv("GO_BF_miRNA_INT_counts.tab",sep="\t")