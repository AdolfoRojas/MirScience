library(MetaVolcanoR)

datasets <- list()
for (file in list.files("./", pattern = "_dif_exp_NS.tsv")){
    print(file)
    name <- gsub("_dif_exp_NS.tsv","",file)
    archivo <- read.delim(file,sep ="\t")
    colnames(archivo)[colnames(archivo) =="Row.names"] <- "ID"
    if ("Gene.symbol" %in% colnames(archivo)){
    archivo$X <- NULL
    nam_best_grid_nosp <- name
    assign(nam_best_grid_nosp, archivo)
    #datasets <- append(datasets,get(nam_best_grid_nosp))
    datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)}
}

meta_degs_rem <- rem_mv(diffexp=datasets,
            pcriteria="adj.P.Val",
            foldchangecol='logFC', 
            genenamecol="Gene.symbol",
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
                   pcriteria="adj.P.Val",
                   foldchangecol='logFC',
                   genenamecol="Gene.symbol",
                   geneidcol="ID",
                   pvalue=0.05,
                   foldchange=0, 
                   metathr=0.01,
                   collaps=FALSE,
                   jobname="MetaVolcano", 
                   outputfolder=".",
                   draw='HTML')

head(meta_degs_vote@metaresult, 3)
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