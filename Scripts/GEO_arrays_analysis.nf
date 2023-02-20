#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.GEO_TABLE = "/media/run-projects/Adolfo/MirScience/Prueba_de_codigo/nash_GEO_samples.csv"
params.Project = ""
params.Organism = ""
params.normalize = false

process Dowload_Datasets { 
    conda "/media/run-projects/software/miniconda3/envs/GEO2R" 
    publishDir "${params.outdir}/Raw_data/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 10
    tag "Download datasets $accession"

    input:
    tuple val(accession),  val(GPL),  val(Bioproject)

    output:  
    tuple val(accession),  val(GPL),  val(Bioproject), path("${accession}.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript 
    library(GEOquery)

    if (file.exists("../../../${params.outdir}Raw_data/${accession}.rds")){
        system("cp ../../../${params.outdir}Raw_data/${accession}.rds .")
    } else { 
        gset <- getGEO("${accession}",GSEMatrix =TRUE, AnnotGPL=TRUE)
        saveRDS(gset, file = "${accession}.rds")
    }
    """   
}
process Analyze { 
    conda "/media/run-projects/software/miniconda3/envs/GEO2R" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 5   
    tag "Analyzing $accession"

    input:
    tuple val(accession),  val(GPL),  val(Bioproject), path(rds_file)
    path(Sample_table)

    output:  
    path("*")

    script:
    """
    #!/usr/bin/env Rscript 
    library(GEOquery)
    library(tidyverse)
    library(limma)
    library(umap)

    GEO_TABLE = read.delim("${Sample_table}",sep=",") 
    def_samples = GEO_TABLE[GEO_TABLE\$Condition != "ND",]

    print("${accession}")
    system(paste0("mkdir -p ","${Bioproject[0]}"))
    Conditions <- c("control","Case")
    df_loop <- def_samples[def_samples\$Accession2 == "${accession}",] 
    GLP <- "${GPL[0]}"
    df_loop <- df_loop[c("Accession","Condition")]
        ################################################################
    # map probe IDs to gene symbols
    gset <- readRDS("${rds_file}")
    
    if (length(gset) > 1) idx <- grep(paste0("GPL",GLP), attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    # make proper column names to match toptable 
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    ################################################################
    gset <- gset[,df_loop\$Accession]
    sml <- df_loop\$Condition
    ################################################################
    # LOG transformation if necessary
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
                (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
    ################################################################
    ################################################################
    # assign samples to groups and set up design matrix
    gs <- factor(sml)
    groups <- make.names(c("Control","Case"))
    levels(gs) <- groups
    ################################################################
    gset\$group <- gs
    design <- model.matrix(~group + 0, gset)
    colnames(design) <- levels(gs)
    fit <- lmFit(gset, design)  # fit linear model
    # set up contrasts of interest and recalculate model coefficients
    cts <- paste(groups[1], groups[2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    # compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000,p.value = 0.05)
    #tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
    write.table(tT, file=paste0("${Bioproject[0]}","_dif_exp.tsv"), row.names=F, sep="\t")
    # Visualize and quality control test results.
    # Build histogram of P-values for all genes. Normal test
    # assumption is that most genes are not differentially expressed.
    pdf(paste0("${Bioproject[0]}","/hist.pdf"))
    tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
    hist(tT2\$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",ylab = "Number of genes", main = "P-adj value distribution")
    dev.off()
    # summarize test results as "up", "down" or "not expressed"
    dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
    # Venn diagram of results
    pdf(paste0("${Bioproject[0]}","/vennDiagram.pdf"))
    vennDiagram(dT, circle.col=palette())
    dev.off()
    # create Q-Q plot for t-statistic
    pdf(paste0("${Bioproject[0]}","/qqt.pdf"))
    t.good <- which(!is.na(fit2\$F)) # filter out bad probes
    qqt(fit2\$t[t.good], fit2\$df.total[t.good], main="Moderated t statistic")
    dev.off()
    # volcano plot (log P-value vs log fold change)
    colnames(fit2) # list contrast names
    ct <- 1        # choose contrast of interest

    tT2\$diffexpressed <- "NS"
    tT2\$diffexpressed[tT2\$logFC > 0 & tT2\$adj.P.Val < 0.05] <- "UP"
    tT2\$diffexpressed[tT2\$logFC < 0 & tT2\$adj.P.Val < 0.05] <- "DOWN"

    ggplot(data = tT2, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + #)
    ggtitle("${accession}")+
    geom_point(size = 1.4) + 
    scale_color_manual("DEGs",values=c("blue", "gray", "red")) +theme_light() + theme(title = element_text(size= 14),plot.title = element_text(hjust = 0.5)) #+
    #geom_vline(xintercept=c(-1, 1), col="red") +
    #geom_hline(yintercept=-log10(0.05), col="red")  

    ggsave(paste0("${Bioproject[0]}","/volcanoplot.pdf"))
    write.table(tT2, file=paste0("${Bioproject[0]}","/","${accession}","_dif_exp_NS.tsv"), row.names=F, sep="\t")

    # MD plot (log fold change vs mean log expression)
    # highlight statistically significant (p-adj < 0.05) probes
    pdf(paste0("${Bioproject[0]}","/plotMD.pdf"))
    plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
    abline(h=0)
    dev.off()
    ################################################################
    # General expression data analysis
    ex <- exprs(gset)
    # box-and-whisker plot  ################################################################################# importante para saber normalizacion
    pdf(paste0("${Bioproject[0]}","/boxplot.pdf"))
    ord <- order(gs)  # order samples by group
    palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
    par(mar=c(7,4,2,1))
    title <- paste("${accession}", "/", annotation(gset), sep ="")
    boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
    legend("topleft", groups, fill=palette(), bty="n")
    dev.off()
    # expression value distribution
    pdf(paste0("${Bioproject[0]}","/plotDensities.pdf"))
    par(mar=c(4,4,2,1))
    title <- paste ("${accession}", "/", annotation(gset), " value distribution", sep ="")
    plotDensities(ex, group=gs, main=title, legend ="topright")
    dev.off()
    # UMAP plot (dimensionality reduction)
    ex <- na.omit(ex) # eliminate rows with NAs
    ex <- ex[!duplicated(ex), ]  # remove duplicates
    ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
        #percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaData = as.data.frame(ump\$layout)
    names(pcaData) = c("PC1","PC2")
    p <- ggplot(pcaData, aes(PC1, PC2, color=gset\$group)) +
        ggtitle(paste("${accession}", sep = "")) +          
        theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold"),aspect.ratio=1) +
        geom_point(size=2) +
        #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() + scale_color_manual("Condition", values=c("#56B4E9", "red")) 
    ggsave(paste0("${Bioproject[0]}","/UMAP.pdf"), plot = p)
    # mean-variance trend, helps to see if precision weights are needed
    pdf(paste0("${Bioproject[0]}","/plotSA.pdf"))
    plotSA(fit2, main=paste0("Mean variance trend, ","${accession}"))
    dev.off()
    system("rm Rplots.pdf")
    """   
}
process Normalize_and_Analyze { 
    conda "/media/run-projects/software/miniconda3/envs/GEO2R" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 5   
    tag "Analyzing $accession"

    input:
    tuple val(accession),  val(GPL),  val(Bioproject), path(rds_file)
    path(Sample_table)

    output:  
    path("*")

    script:
    """
    #!/usr/bin/env Rscript 
    library(GEOquery)
    library(tidyverse)
    library(limma)
    library(umap)

    GEO_TABLE = read.delim("${Sample_table}",sep=",") 
    def_samples = GEO_TABLE[GEO_TABLE\$Condition != "ND",]

    print("${accession}")
    system(paste0("mkdir -p ","${Bioproject[0]}"))
    Conditions <- c("control","Case")
    df_loop <- def_samples[def_samples\$Accession2 == "${accession}",] 
    GLP <- "${GPL[0]}"
    df_loop <- df_loop[c("Accession","Condition")]
        ################################################################
    # map probe IDs to gene symbols
    gset <- readRDS("${rds_file}")
    
    if (length(gset) > 1) idx <- grep(paste0("GPL",GLP), attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    # make proper column names to match toptable 
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    ################################################################
    gset <- gset[,df_loop\$Accession]
    sml <- df_loop\$Condition
    ################################################################
    # LOG transformation if necessary
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
                (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
    ################################################################
    exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
    ################################################################
    # assign samples to groups and set up design matrix
    gs <- factor(sml)
    groups <- make.names(c("Control","Case"))
    levels(gs) <- groups
    ################################################################
    gset\$group <- gs
    design <- model.matrix(~group + 0, gset)
    colnames(design) <- levels(gs)
    fit <- lmFit(gset, design)  # fit linear model
    # set up contrasts of interest and recalculate model coefficients
    cts <- paste(groups[1], groups[2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    # compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000,p.value = 0.05)
    #tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
    write.table(tT, file=paste0("${Bioproject[0]}","_Norm_dif_exp.tsv"), row.names=F, sep="\t")
    # Visualize and quality control test results.
    # Build histogram of P-values for all genes. Normal test
    # assumption is that most genes are not differentially expressed.
    pdf(paste0("${Bioproject[0]}","/Norm_hist.pdf"))
    tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
    hist(tT2\$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",ylab = "Number of genes", main = "P-adj value distribution")
    dev.off()
    # summarize test results as "up", "down" or "not expressed"
    dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
    # Venn diagram of results
    pdf(paste0("${Bioproject[0]}","/Norm_vennDiagram.pdf"))
    vennDiagram(dT, circle.col=palette())
    dev.off()
    # create Q-Q plot for t-statistic
    pdf(paste0("${Bioproject[0]}","/Norm_qqt.pdf"))
    t.good <- which(!is.na(fit2\$F)) # filter out bad probes
    qqt(fit2\$t[t.good], fit2\$df.total[t.good], main="Moderated t statistic")
    dev.off()
    # volcano plot (log P-value vs log fold change)
    colnames(fit2) # list contrast names
    ct <- 1        # choose contrast of interest

    tT2\$diffexpressed <- "NS"
    tT2\$diffexpressed[tT2\$logFC > 0 & tT2\$adj.P.Val < 0.05] <- "UP"
    tT2\$diffexpressed[tT2\$logFC < 0 & tT2\$adj.P.Val < 0.05] <- "DOWN"

    ggplot(data = tT2, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + #)
    ggtitle("${accession}")+
    geom_point(size = 1.4) + 
    scale_color_manual("DEGs",values=c("blue", "gray", "red")) +theme_light() + theme(title = element_text(size= 14),plot.title = element_text(hjust = 0.5)) #+
    #geom_vline(xintercept=c(-1, 1), col="red") +
    #geom_hline(yintercept=-log10(0.05), col="red")  

    ggsave(paste0("${Bioproject[0]}","/Norm_volcanoplot.pdf"))
    write.table(tT2, file=paste0("${Bioproject[0]}","/","Norm_${Bioproject[0]}","_dif_exp_NS.tsv"), row.names=F, sep="\t")

    # MD plot (log fold change vs mean log expression)
    # highlight statistically significant (p-adj < 0.05) probes
    pdf(paste0("${Bioproject[0]}","/Norm_plotMD.pdf"))
    plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
    abline(h=0)
    dev.off()
    ################################################################
    # General expression data analysis
    ex <- exprs(gset)
    # box-and-whisker plot  ################################################################################# importante para saber normalizacion
    pdf(paste0("${Bioproject[0]}","/Norm_boxplot.pdf"))
    ord <- order(gs)  # order samples by group
    palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
    par(mar=c(7,4,2,1))
    title <- paste("${accession}", "/", annotation(gset), sep ="")
    boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
    legend("topleft", groups, fill=palette(), bty="n")
    dev.off()
    # expression value distribution
    pdf(paste0("${Bioproject[0]}","/Norm_plotDensities.pdf"))
    par(mar=c(4,4,2,1))
    title <- paste ("${accession}", "/", annotation(gset), " value distribution", sep ="")
    plotDensities(ex, group=gs, main=title, legend ="topright")
    dev.off()
    # UMAP plot (dimensionality reduction)
    ex <- na.omit(ex) # eliminate rows with NAs
    ex <- ex[!duplicated(ex), ]  # remove duplicates
    ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
        #percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaData = as.data.frame(ump\$layout)
    names(pcaData) = c("PC1","PC2")
    p <- ggplot(pcaData, aes(PC1, PC2, color=gset\$group)) +
        ggtitle(paste("${accession}", sep = "")) +          
        theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold"),aspect.ratio=1) +
        geom_point(size=2) +
        #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() + scale_color_manual("Condition", values=c("#56B4E9", "red")) 
    ggsave(paste0("${Bioproject[0]}","/Norm_UMAP.pdf"), plot = p)
    # mean-variance trend, helps to see if precision weights are needed
    pdf(paste0("${Bioproject[0]}","/Norm_plotSA.pdf"))
    plotSA(fit2, main=paste0("Mean variance trend, ","${accession}"))
    dev.off()
    system("rm Rplots.pdf")
    """   
}
workflow { 
    params.outdir = "Resultados/${params.Organism}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()

    Sample_table = Channel.fromPath(params.GEO_TABLE)
    if (params.normalize == false){    
        print "\nCheck Information:\n $params.Organism\n"
        ST = Channel.fromPath(params.GEO_TABLE).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.taxon == params.Organism)
    } else if (params.normalize == true){
         print "\nCheck Information:\n$params.Project,$params.Organism\n"
         ST = Channel.fromPath(params.GEO_TABLE).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.BioProjectID == params.Project).filter(row -> row.taxon == params.Organism)
         ST.ifEmpty { exit 1, "ERROR: $params.Project not in table or not defined\n" }}
    
    Projects = ST.map { row-> tuple(row.Accession, row.Title, row.Accession2, row.GPL, row.BioProjectID) } 
    Info = Projects.map { it[0,1,2,3] }.unique()
    Runs = Projects.map { it[2,3,4] }.unique().groupTuple()
    Runs.view()
    print 'If something is wrong press "ctrl + C" and change parameters'
    print "Workflow project outdir " + params.outdir
    Dowload_Datasets(Runs)
    if (params.normalize == false){ 
        Analyze(Dowload_Datasets.out,Sample_table)
    } else if (params.normalize == true){
       Normalize_and_Analyze(Dowload_Datasets.out,Sample_table)
    }
}