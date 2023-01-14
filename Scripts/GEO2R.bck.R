################################################################
library(affy)
library(GEOquery)
library(tidyverse)
library(RSQLite)
library(AnnotationDbi)
library(limma)
library(umap)
library(maptools)
GEO_TABLE = read.delim("/media/run-projects/Adolfo/MirScience/Datasets/BioProject/nash_GEO_affi_samples.csv",sep=",") 
def_samples = GEO_TABLE[GEO_TABLE$Condition != "ND",]

for (DataSet in unique(def_samples$Accession2)){ 
  ################################################################
  Conditions <- c("control","Case")
  df_loop <- def_samples[def_samples$Accession2 == DataSet,]
  #DataSet <- "GSE148537" #"GSE200409"  
  GLP <- df_loop$GPL[1]
  df_loop <- df_loop[c("Accession","Condition")]
  ################################################################
  # get supplementary files
  getGEOSuppFiles(DataSet)
  # untar files
  untar(paste0(DataSet,"/",DataSet,"_RAW.tar"), exdir = 'data/')
  # reading in .cel files
  raw.data <- ReadAffy(celfile.path = "data/")
  # performing RMA normalization
  normalized.data <- rma(raw.data) ## ‘RSQLite’, ‘AnnotationDbi’
  # get expression estimates
  ################################################################
  # LOG transformation if necessary
  ex <- exprs(normalized.data)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(normalized.data) <- log2(ex) }
  ################################################################
  # map probe IDs to gene symbols
  gset <- getGEO(DataSet, GSEMatrix = TRUE)
  if (length(gset) > 1) idx <- grep(paste0("GPL",GLP), attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  feature.data <- gset@featureData@data# fetch feature data to get ID - gene symbol mapping
  feature.data <- feature.data[,c(1,11)] # subset
  normalized.expr <- as.data.frame(exprs(normalized.data))
  normalized.expr <- normalized.expr %>%
    rownames_to_column(var = 'ID') %>%
    inner_join(., feature.data, by = 'ID')
  colnames(normalized.expr) <- make.names(colnames(normalized.expr))
  rownames(normalized.expr) <- normalized.expr$ID
  normalized.expr$ID <- NULL
  normalized.expr$Gene.Symbol <- gsub(" /// ","|",normalized.expr$Gene.Symbol)
  normalized.expr <- normalized.expr[c("Gene.Symbol")]
  ################################################################
  colnames(normalized.data) <- gsub("\\.","_",colnames(normalized.data))
  colnames(normalized.data) <- sapply(strsplit(colnames(normalized.data),"_"), "[", 1)
  ################################################################
  normalized.data <- normalized.data[,df_loop$Accession]
  sml <- df_loop$Condition
  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("Control","Case"))
  levels(gs) <- groups
  ################################################################
  normalized.data$group <- gs
  design <- model.matrix(~group + 0, normalized.data)
  colnames(design) <- levels(gs)
  fit <- lmFit(normalized.data, design)  # fit linear model
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000,p.value = 0.05)
  tT <- merge(tT,normalized.expr,by = 0)
  #tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
  write.table(tT, file=paste0(DataSet,"_dif_exp.tsv"), row.names=F, sep="\t")
  # Visualize and quality control test results.
  # Build histogram of P-values for all genes. Normal test
  # assumption is that most genes are not differentially expressed.
  pdf(paste0(DataSet,"/hist.pdf"))
  tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",ylab = "Number of genes", main = "P-adj value distribution")
  dev.off()
  # summarize test results as "up", "down" or "not expressed"
  dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
  # Venn diagram of results
  pdf(paste0(DataSet,"/vennDiagram.pdf"))
  vennDiagram(dT, circle.col=palette())
  dev.off()
  # create Q-Q plot for t-statistic
  pdf(paste0(DataSet,"/qqt.pdf"))
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
  dev.off()
  # volcano plot (log P-value vs log fold change)
  colnames(fit2) # list contrast names
  ct <- 1        # choose contrast of interest
  pdf(paste0(DataSet,"/volcanoplot.pdf"))
  volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
  dev.off()
  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  pdf(paste0(DataSet,"/plotMD.pdf"))
  plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  abline(h=0)
  dev.off()
  ################################################################
  # General expression data analysis
  ex <- exprs(normalized.data)
  # box-and-whisker plot  ################################################################################# importante para saber normalizacion
  pdf(paste0(DataSet,"/boxplot.pdf"))
  ord <- order(gs)  # order samples by group
  palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
  par(mar=c(7,4,2,1))
  title <- paste(DataSet, "/", annotation(normalized.data), sep ="")
  boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
  legend("topleft", groups, fill=palette(), bty="n")
  dev.off()
  # expression value distribution
  pdf(paste0(DataSet,"/plotDensities.pdf"))
  par(mar=c(4,4,2,1))
  title <- paste (DataSet, "/", annotation(normalized.data), " value distribution", sep ="")
  plotDensities(ex, group=gs, main=title, legend ="topright")
  dev.off()
  # UMAP plot (dimensionality reduction)
  pdf(paste0(DataSet,"/UMAP.pdf"))
  ex <- na.omit(ex) # eliminate rows with NAs
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = 2, random_state = 123)
  par(mar=c(3,3,2,6), xpd=TRUE)
  plot(ump$layout, main="UMAP plot, nbrs=2", xlab="", ylab="", col=gs, pch=20, cex=1.5)
  legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
  col=1:nlevels(gs), title="Group", pt.cex=1.5)
  library("maptools")  # point labels without overlaps
  pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
  dev.off()
  # mean-variance trend, helps to see if precision weights are needed
  pdf(paste0(DataSet,"/plotSA.pdf"))
  plotSA(fit2, main=paste0("Mean variance trend, ",DataSet))
  dev.off()
  system("rm -r data")}
  ################
  ##########