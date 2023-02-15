#!/usr/bin/env Rscript
library("CEMiTool")
system("mkdir -p normal_vs_tumoral")
PP_int_df <- read.delim("2_interacciones_sobre_700_gene_symbol_in_GTEx.tab", sep = "\t") #Gene-level_mRNAs_DESeq_norm_matrix_WT.csv
expr0 <- read.delim("mRNAs_TCGA-LUAD_DESeq_norm_matrix_WT.csv", sep = ",")
rownames(expr0) <- sapply(strsplit(as.character(rownames(expr0)),'\\.'), "[", 1)
#expr0 <- read.delim("../analisis_muestras_y_DE/Gene-level_mRNAs_DESeq_norm_matrix_WT.csv", sep = ",")
if (file.exists("2_gene_ID_to_gene_symbol_TCGA.tab") == F) {
        gencode_v38 <- read.delim("gencode_genes.v38.annotation.tab", sep = "\t")
        gencode_v38 <- gencode_v38[c("gene_id","gene_name")]
        gencode_v38 <- gencode_v38[!duplicated(gencode_v38),]
        gencode_v38 <- gencode_v38[gencode_v38$gene_id %in% rownames(expr0), ]
        gencode_v38_sorted <- gencode_v38[order(gencode_v38$gene_id), ]

        for (duplicado in levels(as.factor(gencode_v38_sorted[duplicated(gencode_v38_sorted$gene_name),]$gene_name))){
                count = 0
                for (elementos in rownames(gencode_v38_sorted[gencode_v38_sorted$gene_name == duplicado,])){ 
                        if(count != 0) {
                                gencode_v38_sorted[elementos,]$gene_name <- paste(duplicado, count, sep = "_")
                                print(gencode_v38_sorted[elementos, ]$gene_name)}
                        count = count + 1}}
        write.table(gencode_v38_sorted, file = "2_gene_ID_to_gene_symbol_TCGA.tab", row.names = F, quote = F, col.names = T, sep = "\t")}

gencode_v38_sorted <- read.delim("2_gene_ID_to_gene_symbol_TCGA.tab", sep = "\t")
for (index in rownames(expr0)){
    print(index)
    if (any(gencode_v38_sorted$gene_id == index)){
        rownames(expr0)[rownames(expr0) == index] <- gencode_v38_sorted[gencode_v38_sorted$gene_id == index, ]$gene_name}}

PP_int_df <- merge(PP_int_df, gencode_v38_sorted, by.x = "protein1", by.y = "gene_id")
PP_int_df <- PP_int_df[c("gene_name","protein2")]
colnames(PP_int_df)[colnames(PP_int_df) == "gene_name"] <- "protein1"
PP_int_df <- merge(PP_int_df, gencode_v38_sorted, by.x = "protein2", by.y = "gene_id")
PP_int_df <- PP_int_df[c("gene_name","protein1")]
colnames(PP_int_df)[colnames(PP_int_df) == "gene_name"] <- "protein2"
colnames(expr0) <- gsub("\\.", "-", colnames(expr0))



sample_annot <- read.delim("sample_mRNA_info_TCGA-LUAD.tab", sep = "\t")
gmt_info <- read_gmt("h.all.v7.4.symbols.gmt")
expr0 <- expr0[sample_annot$column]
sample_annot$Class <- sample_annot$Condition
sample_annot$SampleName <- sample_annot$column
#sample_annot$Class[sample_annot$Class == "NT"] <- "Normal control"
#sample_annot$Class[sample_annot$Class == "TP"] <- "Tumoral"
sample_annot <- sample_annot[c("SampleName","Class")]

if (identical(colnames(expr0), sample_annot$SampleName) == T){
        cem <- cemitool(expr0, sample_annot, interactions=PP_int_df, gmt_info, filter= T, plot = TRUE, verbose=T, apply_vst = T, gsea_max_size=3000) # 
        
        write_files(cem, directory="normal_vs_tumoral/Tables", force = T) # write analysis results into files
        save_plots(cem, "all", directory="normal_vs_tumoral/Plots", force = T)# save all plots
} else {
        print("problema de ejecucion")}


library(clusterProfiler)
module_interactions_info <- read.delim("normal_vs_tumoral/Tables/module.tsv", sep = "\t")
module_interactions_info <- read.delim("DE_analisis/mRNAs_TCGA-LUAD_Normal_vs_Tumoral_DE.csv", sep = ",")[c("gene_name","log2FoldChange")]
module_interactions_info$diffexpressed <- "NS"
module_interactions_info$diffexpressed[module_interactions_info$log2FoldChange > 2] <- "UP"
module_interactions_info$diffexpressed[module_interactions_info$log2FoldChange < -2] <- "DOWN"
module_interactions_info$log2FoldChange <- NULL
colnames(module_interactions_info) <- c("genes","modules")
write.table(module_interactions_info[module_interactions_info$modules != "NS",]$genes, sep = "\t", file = "DEGs_list.txt", row.names = F, quote = F, col.names = F) 
write.table(module_interactions_info[module_interactions_info$modules == "DOWN",]$genes, sep = "\t", file = "DEGs_DOWN_list.txt", row.names = F, quote = F, col.names = F) 
write.table(module_interactions_info[module_interactions_info$modules == "UP",]$genes, sep = "\t", file = "DEGs_UP_list.txt", row.names = F, quote = F, col.names = F) 
module_interactions_info <- module_interactions_info[module_interactions_info$modules != "NS",]


if (file.exists("MSigDB_Hallmark_2020.gmt") == F) {
        url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
        download.file(url, destfile = "MSigDB_Hallmark_2020.gmt")}

MSigDB_Hallmark_2020 <- read.gmt("MSigDB_Hallmark_2020.gmt")
Comparacion <- compareCluster(genes~modules, data=module_interactions_info,fun=enricher, TERM2GENE=MSigDB_Hallmark_2020)
png(paste("MSigDB_Hallmark_2020","_Module_enrichment.png",sep=""),width = 720, height=720)    
dotplot(Comparacion)
dev.off()
system("scp -P 1313 MSigDB_Hallmark_2020_Module_enrichment.png adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/TCGA/Lung/normal_vs_tumoral/")

if (file.exists("KEGG_2021_Human.gmt") == F) {
        url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"
        download.file(url, destfile = "KEGG_2021_Human.gmt")}

KEGG_2021_Human <- read.gmt("KEGG_2021_Human.gmt")
Comparacion <- compareCluster(genes~modules, data=module_interactions_info,fun=enricher, TERM2GENE=KEGG_2021_Human)
png(paste("KEGG_2021_Human","_Module_enrichment.png",sep=""),width = 720, height=720)    
dotplot(Comparacion)
dev.off()
system("scp -P 1313 KEGG_2021_Human_Module_enrichment.png adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/TCGA/Lung/normal_vs_tumoral/")

if (file.exists("WikiPathway_2021_Human.gmt") == F) {
        url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathway_2021_Human"
        download.file(url, destfile = "WikiPathway_2021_Human.gmt")}

WikiPathway_2021_Human <- read.gmt("WikiPathway_2021_Human.gmt")
Comparacion <- compareCluster(genes~modules, data=module_interactions_info,fun=enricher, TERM2GENE=WikiPathway_2021_Human)
png(paste("WikiPathway_2021_Human","_Module_enrichment.png",sep=""),width = 720, height=720)    
dotplot(Comparacion)
dev.off()
system("scp -P 1313 WikiPathway_2021_Human_Module_enrichment.png adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/TCGA/Lung/normal_vs_tumoral/")

