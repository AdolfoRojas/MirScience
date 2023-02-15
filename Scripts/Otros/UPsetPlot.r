#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)
library(UpSetR)
colors <- clusterExperiment::bigPalette
datasets <- list()
for (file in list.files("./", pattern = "_DE_NS.csv")){
    print(file)
    name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
    archivo <- read.delim(file,sep =",")
    archivo <- archivo[archivo$padj <= 0.05,]
    colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
    datasets[[name]] <- archivo$gene_id}

png(file="filename.png", res=400, width= 16, height= 9,units = "in") # or other device
upset(fromList(datasets), order.by = "freq",nintersects = 70, nsets = length(datasets))
dev.off()

seurat_theme <- function(){ # Just aesthetics configuration of plots, defined by carlos
  theme_bw() +
    theme(panel.background = element_rect(colour = "black", size=0.1), 
          plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=14),legend.text = element_text(size = 14), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())}

Datos<- NULL
contador <-0
for (file in list.files("./", pattern = "_DE_NS.csv")){
        print(file)
        name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
        archivo <- read.delim(file,sep=",")
        archivo$Project <- name
        archivo$X <- NULL
        archivo <- archivo[!is.na(archivo$padj),]
        archivo <- archivo[archivo$padj <= 0.05,]
        colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
        if (contador==0){
            Datos <- archivo
        } else{
            Datos<- rbind(Datos,archivo)
        }
        contador <-contador+1
}

Datos$important_labels <- ""
for (i in unique(Datos$Project)){
  Datos_in_loop <- Datos[Datos$Project == i,]
  important_genes_DOWN <- head(Datos_in_loop[order(Datos_in_loop$log2FoldChange),]$gene_name,5)
  important_genes_UP <- tail(Datos_in_loop[order(Datos_in_loop$log2FoldChange),]$gene_name,5)
  for (k in c(important_genes_DOWN,important_genes_UP)){
    Datos[Datos$Project == i & Datos$gene_name == k,]$important_labels <- k
  }
}

pos <- position_jitter(width = 0.3, seed = 2)
ggplot(Datos,aes(x = Project, y = log2FoldChange, fill = Project)) + geom_jitter(aes(color = Project, size = -log(padj)),position = pos) + geom_text_repel(aes(label = Datos$important_labels),box.padding = 0.5,position = pos,max.overlaps = Inf) + theme(axis.text.x=element_blank()) + labs(x = "")+ seurat_theme() + ggtitle(paste0("Projects\n significant DEGs")) + scale_color_manual(values=colors)  + scale_y_continuous(limits = c(-1 * max(round(abs(Datos$log2FoldChange)))-1, max(round(abs(Datos$log2FoldChange))))+1) + geom_hline(yintercept=0) + geom_hline(yintercept=0.25, linetype="dashed") + geom_hline(yintercept=-0.25, linetype="dashed") + theme(axis.text.x=element_blank())
ggsave("Jitterplot_DEGs.png", units = "in",dpi = 300, width = 21, height = 9,limitsize = FALSE, bg = "white")

#####################################################################################################################################################

library(ggplot2)
library(ggrepel)
library(UpSetR)
colors <- clusterExperiment::bigPalette
datasets <- list()

all <- read.delim("All/All_important_genes_votecounting_metaanalisis.tsv",sep ="\t")
mouse <- read.delim("Mus_musculus/Mus_musculus_important_genes_votecounting_metaanalisis.tsv",sep ="\t")
human <- read.delim("Homo_sapiens/Homo_sapiens_important_genes_votecounting_metaanalisis.tsv",sep ="\t")
diopt <- read.delim("/media/run-projects/Adolfo/MirScience/RNAseq/Mus_musculus/diopt8.5_results_Mm27_to_Hs.high.tab", sep= "\t")[c("Search.Term","Ensmbl.ID...link.HPA.")]

datasets[["Inter-espcie"]]  <- all$gene_id
datasets[["Homo_sapiens"]]  <- human$gene_id
datasets[["Mus_musculus"]] <- mouse$gene_id
datasets[["Mus_musculus_orthologs"]] <- mouse$gene_id[mouse$gene_id %in% diopt$Search.Term]
datasets[["Human_orthologs"]] <- human$gene_id[human$gene_id %in% diopt$Ensmbl.ID...link.HPA.]

deg_orthologs <- diopt[diopt$Search.Term %in% mouse$gene_id,]
deg_orthologs_df <- mouse[mouse$gene_id %in% deg_orthologs$Search.Term,]
dim(deg_orthologs_df[deg_orthologs_df$degvcount == "0.Down-regulated",])
dim(deg_orthologs_df[deg_orthologs_df$degvcount == "2.Up-regulated",])

write.table(deg_orthologs$Ensmbl.ID...link.HPA., file="mouse_orth_DEGs.tsv", row.names=F, sep="\t",quote=F)


png(file="Upset_shared_orthologs.png", res=400, width= 21, height= 9,units = "in") # or other device
upset(fromList(datasets), order.by = "freq",nintersects = 70, nsets = length(datasets),set_size.show = TRUE)
dev.off()


python3
from bioinfokit import analys, visuz
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
v = venn2( (27, 26, 119), alpha = 1 , set_labels = ('Human', 'Mouse', 'Orthologs'),set_colors=('#a1b8e0', '#a1b8e0'))
plt.savefig('Venn_orthologs_DEGs2.png', dpi = 500)