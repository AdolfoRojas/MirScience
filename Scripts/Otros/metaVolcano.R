library(MetaVolcanoR)

datasets <- list()

for (file in list.files("Homo_sapiens/", pattern = "_Control_vs_Case_DE_NS.csv")){
    print(file)
    name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
    archivo <- read.delim(paste0("Homo_sapiens/",file),sep =",")
    colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
    archivo$X <- NULL
    nam_best_grid_nosp <- name
    assign(nam_best_grid_nosp, archivo)
    datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)
}
IDs_Hs <- read.delim("IDs_equivalencias_gencode.v38.txt", sep= "\t")
IDs_Mm <- read.delim("IDs_equivalencias_gencode.vM27.txt", sep= "\t")
IDs_Mm$gene_id <- sapply(strsplit(as.character(IDs_Mm$gene_id),'\\.'), "[", 1)
diopt <- read.delim("Mus_musculus/diopt8.5_results_Mm27_to_Hs.high.tab", sep= "\t")[c("Search.Term","Ensmbl.ID...link.HPA.")]
names(diopt) <- c("GeneID_Mm","GeneID_Hs")
diopt <- merge(diopt,IDs_Hs,by.x="GeneID_Hs", by.y="gene_id")

for (file in list.files("Mus_musculus/", pattern = "_Control_vs_Case_DE_NS.csv")){
    print(file)
    name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
    archivo <- read.delim(paste0("Mus_musculus/",file),sep =",")
    colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
    archivo$gene_name <- NULL
    archivo <- merge(archivo,diopt,by.x="gene_id",by.y="GeneID_Mm")
    archivo$gene_id <- NULL
    names(archivo)[names(archivo)=="GeneID_Hs"] <- "gene_id"
    archivo$X <- NULL
    archivo <- archivo[order(archivo$padj),]
    archivo <- archivo[!duplicated(archivo$gene_id),]
    nam_best_grid_nosp <- name
    assign(nam_best_grid_nosp, archivo)
    datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)
}


for (file in list.files("./", pattern = "_Control_vs_Case_DE_NS.csv")){
    print(file)
    name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
    archivo <- read.delim(file,sep =",")
    colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
    archivo$X <- NULL
    nam_best_grid_nosp <- name
    assign(nam_best_grid_nosp, archivo)
    #datasets <- append(datasets,get(nam_best_grid_nosp))
    datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)
}





meta_degs_rem <- rem_mv(diffexp=datasets,
            pcriteria="padj",
            foldchangecol='log2FoldChange', 
            genenamecol='gene_name',
            #geneidcol="gene_id",
            collaps=TRUE,
            #llcol='CI.L',
            #rlcol='CI.R',
            vcol="lfcSE", 
            cvar=FALSE,
            metathr=0.01,
            jobname="MetaVolcano",
            outputfolder=".", 
            draw='HTML',
            ncores=30)

head(meta_degs_rem@metaresult, 3)
meta_degs_rem@MetaVolcano
#draw_forest(remres=meta_degs_rem,
#        gene="STMN2",
#        genecol="gene_name",
#        foldchangecol="log2FoldChange",
#        #llcol="randomCi.lb", 
#        #rlcol="randomCi.ub",
#        jobname="MetaVolcano",
#        outputfolder=".",
#        draw="PDF")

#draw_forest(remres=meta_degs_rem,
#        gene="COL6A6",
#        genecol="Symbol", 
#        foldchangecol="Log2FC",
#        llcol="CI.L", 
#        rlcol="CI.R",
#        jobname="MetaVolcano",
#        outputfolder=".",
#        draw="PDF")

###############################################

meta_degs_vote <- votecount_mv(diffexp=datasets,
                   pcriteria='padj',
                   foldchangecol='log2FoldChange',
                   genenamecol='gene_name',
                   geneidcol="gene_id",
                   pvalue=0.05,
                   foldchange=0, 
                   metathr=0.01,
                   collaps=FALSE,
                   jobname="MetaVolcano", 
                   outputfolder=".",
                   draw='HTML')

archivo[duplicated(archivo$gene_id),]
archivo[order(archivo$padj),]

head(meta_degs_vote@metaresult, 3)

#important_genes <- meta_degs_vote@metaresult[meta_degs_vote@metaresult$ndeg >= (length(datasets)-3) & meta_degs_vote@metaresult$ndeg == abs(meta_degs_vote@metaresult$ddeg) & meta_degs_vote@metaresult$ndeg > 1,]
important_genes <- meta_degs_vote@metaresult[meta_degs_vote@metaresult$degvcount != "1.Unperturbed",]


library(ggplot2)
# Basic histogram
ggplot(important_genes, aes(x=ndeg)) + geom_histogram()
ggsave("histograma.pdf")

pdf("histograma.pdf")
hist(important_genes$ndeg)
dev.off()

if (grepl("ENSMUSG",important_genes$gene_id[1])){
        important_genes <- merge(IDs_Mm,important_genes,by="gene_id")
} else {
        important_genes <- merge(IDs_Hs,important_genes,by="gene_id")
}



write.table(important_genes, file="important_genes_votecounting_metaanalisis.tsv", row.names=F, sep="\t")
write.table(important_genes$gene_id, file="important_genes_votecounting_metaanalisis_GeneID.txt", row.names=F, sep="\t",quote=F)
#write.table(IDs_Mm$gene_id, file="GeneIDs_Mm.txt", row.names=F, sep="\t",quote=F)
meta_degs_vote@degfreq

meta_degs_vote@MetaVolcano

###############################################

meta_degs_comb <- combining_mv(diffexp=datasets,
                   pcriteria='padj', 
                   foldchangecol='log2FoldChange',
                   genenamecol='gene_name',
                   #geneidcol="gene_id",
                   metafc='Mean',
                   metathr=0.01, 
                   collaps=TRUE,
                   jobname="MetaVolcano",
                   outputfolder=".",
                   draw='HTML')

head(meta_degs_comb@metaresult, 3)
meta_degs_comb@MetaVolcano




    IDs_Hs <- read.delim("../IDs_equivalencias_gencode.v38.txt", sep= "\t")
    IDs_Mm <- read.delim("../IDs_equivalencias_gencode.vM27.txt", sep= "\t")
    diopt <- read.delim("diopt8.5_results_Mm27_to_Hs.high.tab", sep= "\t")

    IDs_Mm\$gene_id <- sapply(strsplit(as.character(IDs_Mm\$gene_id),'\\.'), "[", 1)

    if (grepl("ENSMUSG",important_genes\$gene_id[1])){
            important_genes <- merge(IDs_Mm,important_genes,by="gene_id")
    } else {
            important_genes <- merge(IDs_Hs,important_genes,by="gene_id")
    }