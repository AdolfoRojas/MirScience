library(pdftools)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
library(stringr)
register(MulticoreParam(15))
expr0 <- read.delim("Matrix_count.txt", sep = "\t")
colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
colnames(expr0) <- gsub("-sorted-bam", "", colnames(expr0))
rownames(expr0) <- sapply(strsplit(expr0$Geneid,"\\."), "[", 1)
expr0$Geneid <- NULL
expr0$Chr <- NULL
expr0$Length <- NULL
expr0$Start <- NULL
expr0$End <- NULL
expr0$Strand <- NULL
gencode_info <- read.delim("IDs_equivalencias_gencode.vM27.txt", sep = "\t")
gencode_info$gene_id <- sapply(strsplit(gencode_info$gene_id,"\\."), "[", 1)

sample_data <- read.delim("PRJNA763814.csv", sep = ",")
sample_data$Title_y <- NULL
names(sample_data) <- c("column", "Condition")
coldata <- sample_data 
matrix <- as.data.frame(expr0)    
RNAs <- "mRNAs"    
cts <- matrix[coldata$column]     
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
dds <- dds[keep,]
dds$Condition <- factor(dds$Condition, levels = c("Control","Case"))
dds$Condition <- droplevels(dds$Condition)
dds <- DESeq(dds, fitType="local", parallel = TRUE)        
DESeq_norm <- counts(dds, normalized=T)
write.table(DESeq_norm, file = paste0(RNAs,"_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",")  
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  ggtitle(paste(RNAs, sep = "")) +          
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("#56B4E9", "red"))
ggsave(paste(RNAs,"_PCA_plot.tiff", sep = ""), plot = p, width = 8, height = 8, dpi = 300, units = "in")
res <- results(dds, parallel = TRUE)
res2 <- as.data.frame(res)
if (RNAs == "mRNAs") {
res2 <- merge(res2, gencode_info, by.x = "row.names",by.y = "gene_id")
}
res2$L2FC_ABS <- abs(res2$log2FoldChange)
res2 <- res2[order(res2$L2FC_ABS,decreasing=T),]
write.csv(na.omit(res2[res2$padj < 0.05,]), file=paste(RNAs,"_",levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.csv", sep= ""))
write.csv(res2, file=paste(RNAs,"_",levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE_NS.csv", sep= ""))
        