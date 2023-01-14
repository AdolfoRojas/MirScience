#!/usr/bin/env Rscript
# Authors: Adolfo Rojas, Michel Naslavsky, Vinicius Maracaja
# Affiliation: LIB, Facultad de Ciencias Quimicas y Farmaceuticas, Universidad de Chile, Universidad de Sao Paulo.
library(pdftools)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
library(TCGAbiolinks)
library(SummarizedExperiment)
#-----------------------------------------------------------------------------------------------------------------------------------------#
Cancer <- "LUAD"            ##Enter Cancer study abbreviation LUAD PAAD

GDCquery_clinic(paste("TCGA-", Cancer, sep = ""), type = "clinical", save.csv = T) # Descarga las caracteristicas clinicas de las muestras
Clinical_data <- read.table(paste("TCGA-",  Cancer, "_clinical.csv", sep = ""), header=T, sep = ",")
Reduced_Clinical_data <- Clinical_data[c("submitter_id", "primary_diagnosis", "ethnicity", "race", "gender","age_at_diagnosis")] # Caracteristicas Clinicaas de Interes
length(unique(Reduced_Clinical_data$submitter_id))

CancerProject <- paste0("TCGA-",Cancer)
DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","mRNA_gene_quantification",".rda")
query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", legacy = FALSE) #HTSeq - Counts
samplesDown.miR <- getResults(query.miR,cols=c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
tumor_sample_df <- as.data.frame(dataSmTP.miR)
colnames(tumor_sample_df) <- "Sample"
tumor_sample_df$Condition <- "Tumoral" 
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
normal_sample_df <- as.data.frame(dataSmNT.miR)
colnames(normal_sample_df) <- "Sample"
normal_sample_df$Condition <- "Normal" 
normal_sample_df$Patient <- substring(normal_sample_df$Sample,1,12)
tumor_sample_df$Patient <- substring(tumor_sample_df$Sample,1,12)
normal_sample_df$Sample <- substring(normal_sample_df$Sample,1,15)
tumor_sample_df$Sample <- substring(tumor_sample_df$Sample,1,15)
mRNA_df_info <- merge(normal_sample_df,tumor_sample_df, by= "Patient")

queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
GDCdownload(query = queryDown.miR, directory = DataDirectory)
dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)
matrix <- assay(dataAssy.miR,"unstranded")
write.table(matrix, sep = "\t", file = paste0("mRNAs_gene-level_Counts_",CancerProject,"_matrix.tab"), row.names = T, quote = F, col.names = T)
#######################################################################################################################################
#              miRNAs  contiene el script tesis/1_expression_data/TCGA_data/5_Get_miRNAs_count_matrix.R
#######################################################################################################################################
FileNameData <- paste0(DataDirectory, "_","miRNA_gene_quantification",".rda")
query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE)
samplesDown.miR <- getResults(query.miR,cols=c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
tumor_sample_df <- as.data.frame(dataSmTP.miR)
colnames(tumor_sample_df) <- "Sample"
tumor_sample_df$Condition <- "Tumoral" 
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
normal_sample_df <- as.data.frame(dataSmNT.miR)
colnames(normal_sample_df) <- "Sample"
normal_sample_df$Condition <- "Normal" 
normal_sample_df$Patient <- substring(normal_sample_df$Sample,1,12)
tumor_sample_df$Patient <- substring(tumor_sample_df$Sample,1,12)
normal_sample_df$Sample <- substring(normal_sample_df$Sample,1,15)
tumor_sample_df$Sample <- substring(tumor_sample_df$Sample,1,15)
miRNA_df_info <- merge(normal_sample_df,tumor_sample_df, by= "Patient")

queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
GDCdownload(query = queryDown.miR, directory = DataDirectory)
dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)
rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID
read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))] # using read_count's data
matrix <- dataAssy.miR[, read_countData]
colnames(matrix) <- gsub("read_count_","", colnames(matrix))
names <- rownames(dataAssy.miR)
matrix<- cbind(names,matrix)      
write.table(matrix, sep = "\t", file = paste0("miRNAs_counts_",CancerProject,".tab"), row.names = F, quote = F, col.names = T) 
#######################################################################################################################################
#miRNA_df_info = miRNA_df_info[miRNA_df_info$Patient %in% mRNA_df_info$Patient,]
#mRNA_df_info = mRNA_df_info[mRNA_df_info$Patient %in% miRNA_df_info$Patient,]
miRNA_df_info = miRNA_df_info[!duplicated(miRNA_df_info$Patient),]
mRNA_df_info = mRNA_df_info[!duplicated(mRNA_df_info$Patient),]

sample_mRNA_info1 = mRNA_df_info[c("Patient","Sample.x","Condition.x")]
colnames(sample_mRNA_info1) = c("Patient","column","Condition")
sample_mRNA_info2 = mRNA_df_info[c("Patient","Sample.y","Condition.y")]
colnames(sample_mRNA_info2) = c("Patient","column","Condition")
sample_mRNA_info = rbind(sample_mRNA_info1, sample_mRNA_info2)

sample_miRNA_info1 = miRNA_df_info[c("Patient","Sample.x","Condition.x")]
colnames(sample_miRNA_info1) = c("Patient","column","Condition")
sample_miRNA_info2 = miRNA_df_info[c("Patient","Sample.y","Condition.y")]
colnames(sample_miRNA_info2) = c("Patient","column","Condition")
sample_miRNA_info = rbind(sample_miRNA_info1, sample_miRNA_info2)
write.table(sample_miRNA_info, sep = "\t", file = paste0("sample_miRNA_info_",CancerProject,".tab"), row.names = F, quote = F, col.names = T) 
write.table(sample_mRNA_info, sep = "\t", file = paste0("sample_mRNA_info_",CancerProject,".tab"), row.names = F, quote = F, col.names = T) 
# sample_miRNA_info = rbind(normal_sample_df,tumor_sample_df)[c("Sample","Condition")]
# colnames(sample_miRNA_info) = c("column","Condition")

system("mkdir -p DE_analisis")
for (RNAs in c("mRNAs","miRNAs")) {
        print(RNAs)        
        if (RNAs == "mRNAs") {
          expr0 <- read.delim(paste0("mRNAs_gene-level_Counts_",CancerProject,"_matrix.tab"), sep = "\t")
          coldata <- sample_mRNA_info[c("column","Condition")]           
        } else {
          expr0 <- read.delim(paste0("miRNAs_counts_",CancerProject,".tab"), sep = "\t", row.names=1)  # expr0 <- read.delim("matrix_miRNAs_isoforms_with_names.tab", sep = "\t", row.names=2)
          coldata <- sample_miRNA_info[c("column","Condition")]} # sample_miRNA_info <- read.delim("IDs_TCGA_samples_miRNAs_isoforms.tab",sep = "\t")
        #sample_miRNA_info <- sample_miRNA_info[c("cases","sample_type")]        
        #colnames(sample_miRNA_info) <- c("column","Condition")
        #sample_miRNA_info[sample_miRNA_info$Condition == "Primary Tumor",]$Condition <- "Tumoral"
        #sample_miRNA_info[sample_miRNA_info$Condition == "Solid Tissue Normal",]$Condition <- "Normal"
        #sample_miRNA_info <- sample_miRNA_info[sample_miRNA_info$Condition == "Tumoral" | sample_miRNA_info$Condition == "Normal",]
        #sample_miRNA_info$Patient <- substring(sample_miRNA_info$column, 1, 12) 
        #tumor <- sample_miRNA_info[sample_miRNA_info$Condition=="Tumoral",]
        #normal <- sample_miRNA_info[sample_miRNA_info$Condition=="Normal",]
        #normal <- normal[!duplicated(normal$Patient),]
        #tumor <- tumor[!duplicated(tumor$Patient),]
        #tumor <- tumor[tumor$Patient %in% normal$Patient,]
        #sample_miRNA_info <- rbind(normal,tumor)       
        # sample_miRNA_info$column <- substring(sample_miRNA_info$column, 1, 15) 
        #coldata <- sample_miRNA_info[c("column","Condition")]
        rownames(coldata) <- coldata$column
        colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
        colnames(expr0) <- substring(colnames(expr0), 1, 15)        
        matrix <- as.data.frame(expr0)       
        cts <- matrix[coldata$column]     
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(cts), colData = coldata, design = ~ Condition)
        keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
        dds <- dds[keep,]
        dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
        dds$Condition <- droplevels(dds$Condition)
        dds <- DESeq(dds, fitType="local", parallel = TRUE)        
        DESeq_norm <- counts(dds, normalized=T)
        write.table(DESeq_norm, file = paste0(RNAs,"_",CancerProject,"_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",") # write.table(DESeq_norm, file = "miRNAs_isoform_lung_DESeq_norm_matrix_WT.csv", row.names = T, quote = F, col.names = T, sep = ",")
        vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
        pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
          ggtitle(paste(CancerProject, "\n",RNAs, sep = "")) +          
          theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
          geom_point(size=2) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          coord_fixed() +
          scale_color_manual(values=c("#56B4E9", "red"))
        ggsave(paste("DE_analisis/",RNAs,"_",CancerProject,"_PCA_plot.tiff", sep = ""), plot = p, width = 8, height = 8, dpi = 300, units = "in")
        res <- results(dds, parallel = TRUE)
        res2 <- as.data.frame(res)
        if (RNAs == "mRNAs") {
          gencode <- read.delim("/media/storage/Adolfo/MirScience/gencode_genes.v38.annotation.tab",sep= "\t")
          gencode <- gencode[c("gene_id","gene_name")]
          res2$gene_id <- rownames(res2)
          res2$gene_id <- sapply(strsplit(res2$gene_id,"\\."), "[", 1)
          res2 <- merge(res2, gencode, by = "gene_id")
        }
        res2$L2FC_ABS <- abs(res2$log2FoldChange)
        res2 <- res2[order(res2$L2FC_ABS,decreasing=T),]
        write.csv(na.omit(res2[res2$padj < 0.05,]), file=paste("DE_analisis/", RNAs,"_",CancerProject, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.csv", sep= ""))
        write.csv(res2, file=paste("DE_analisis/", RNAs,"_",CancerProject, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE_NS.csv", sep= ""))
        }