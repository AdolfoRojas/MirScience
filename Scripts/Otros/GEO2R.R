################################################################
library(GEOquery)
library(tidyverse)
library(limma)
library(umap)
GEO_TABLE = read.delim("/media/run-projects/Adolfo/MirScience/Datasets/BioProject/nash_GEO_samples.csv",sep=",") 
def_samples = GEO_TABLE[GEO_TABLE$Condition != "ND",]

for (DataSet in unique(def_samples$Accession2)){ 
  ################################################################
  print(DataSet)
  system(paste0("mkdir -p ",DataSet))
  Conditions <- c("control","Case")
  df_loop <- def_samples[def_samples$Accession2 == DataSet,] 
  GLP <- df_loop$GPL[1]
  df_loop <- df_loop[c("Accession","Condition")]
    ################################################################
  # map probe IDs to gene symbols
  gset <- getGEO(DataSet,GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset) > 1) idx <- grep(paste0("GPL",GLP), attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  ################################################################
  gset <- gset[,df_loop$Accession]
  sml <- df_loop$Condition
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
  gset$group <- gs
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

tT2$diffexpressed <- "NS"
tT2$diffexpressed[tT2$logFC > 0 & tT2$adj.P.Val < 0.05] <- "UP"
tT2$diffexpressed[tT2$logFC < 0 & tT2$adj.P.Val < 0.05] <- "DOWN"

ggplot(data = tT2, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + #)
  ggtitle(DataSet)+
  geom_point(size = 1.4) + 
  scale_color_manual("DEGs",values=c("blue", "gray", "red")) +theme_light() + theme(title = element_text(size= 14),plot.title = element_text(hjust = 0.5)) #+
  #geom_vline(xintercept=c(-1, 1), col="red") +
  #geom_hline(yintercept=-log10(0.05), col="red")  

ggsave(paste0(DataSet,"/volcanoplot.pdf"))
write.table(tT2, file=paste0(DataSet,"/",DataSet,"_dif_exp_NS.tsv"), row.names=F, sep="\t")

  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  pdf(paste0(DataSet,"/plotMD.pdf"))
  plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  abline(h=0)
  dev.off()
  ################################################################
  # General expression data analysis
  ex <- exprs(gset)
  # box-and-whisker plot  ################################################################################# importante para saber normalizacion
  pdf(paste0(DataSet,"/boxplot.pdf"))
  ord <- order(gs)  # order samples by group
  palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
  par(mar=c(7,4,2,1))
  title <- paste(DataSet, "/", annotation(gset), sep ="")
  boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
  legend("topleft", groups, fill=palette(), bty="n")
  dev.off()
  # expression value distribution
  pdf(paste0(DataSet,"/plotDensities.pdf"))
  par(mar=c(4,4,2,1))
  title <- paste (DataSet, "/", annotation(gset), " value distribution", sep ="")
  plotDensities(ex, group=gs, main=title, legend ="topright")
  dev.off()
  # UMAP plot (dimensionality reduction)
  ex <- na.omit(ex) # eliminate rows with NAs
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
    #percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData = as.data.frame(ump$layout)
  names(pcaData) = c("PC1","PC2")
  p <- ggplot(pcaData, aes(PC1, PC2, color=gset$group)) +
    ggtitle(paste(DataSet, sep = "")) +          
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold"),aspect.ratio=1) +
    geom_point(size=2) +
    #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + scale_color_manual("Condition", values=c("#56B4E9", "red")) 
  ggsave(paste0(DataSet,"/UMAP.pdf"), plot = p)
  # mean-variance trend, helps to see if precision weights are needed
  pdf(paste0(DataSet,"/plotSA.pdf"))
  plotSA(fit2, main=paste0("Mean variance trend, ",DataSet))
  dev.off()}
  #################################################################################