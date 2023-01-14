#!/usr/bin/env Rscript
library(ggplot2)
library(stringr)
library(tidyr)
library(jpeg)
library(ggpubr)
file <- read.delim("normal_vs_tumoral/Tables/enrichment_nes.tsv")
colnames(file)[1] <- "Modules"
df <- pivot_longer(file, cols=2:3, names_to = "Class", values_to = "NES")
df <- df[str_order(df$Modules, numeric = TRUE),]
df$Modules <- factor(df$Modules,levels=rev(unique(df$Modules)))


ggplot(df, aes(Class, Modules)) +  
  geom_point(aes(size = abs(NES), colour=NES)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() + theme(plot.title = element_text(face = "bold.italic", hjust= -0.1, size = 24)) + 
  labs(size = "NES")

ggsave("Figura_GSEA_co-expresion.png", dpi = 600, units = "in")
system("scp -P 1313 Figura_GSEA_co-expresion.png adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/TCGA/Lung/normal_vs_tumoral/")