#!/usr/bin/env Rscript
# Script generado durante la Tesis
require(TCGAbiolinks)
library(pdftools)
library(DESeq2)
library(ggplot2)
library(mdp)
library(BiocParallel)
register(MulticoreParam(20))
CancerProject <- "TCGA-LUAD"
DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","miRNA_isoform_quantification",".rda")
query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Isoform Expression Quantification", legacy = FALSE)
a <- as.data.frame(query.miR[["results"]][[1]])
a <- a[c("sample_type","file_name","cases")]
write.table(a, file = paste(DataDirectory, "/IDs_TCGA_samples.tab", sep = ""), sep ='\t', row.names = F, quote = F)
samplesDown.miR <- getResults(query.miR,cols=c("cases"))
dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Isoform Expression Quantification", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
GDCdownload(query = queryDown.miR, directory = DataDirectory)
dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)




#!/usr/bin/env Rscript
library(miRBaseConverter)
matrix <- read.delim("matrix_miRNAs_isoforms.tab", sep ='\t')
Accessions=matrix$X
result2=miRNA_AccessionToName(Accessions,targetVersion="v22")
result2 = result2[complete.cases(result2), ]
names(matrix)[1] <- "Accession"
matrix2 <- merge(result2, matrix, by = "Accession")
length(matrix$Accession)
length(matrix2$TargetName)
length(result2$TargetName)
names(matrix2) <- gsub("\\.", "-", names(matrix2))
write.table(matrix2, file = "matrix_miRNAs_isoforms_with_names.tab", sep ='\t', row.names = F, quote = F)