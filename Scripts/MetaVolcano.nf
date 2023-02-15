#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.Dir = "/media/run-projects/Adolfo/MirScience/RNAseq/"
params.Organism = ""
params.Term = ""
params.Combining_approach = false
params.Random_Effect_Model = false
params.Inter_species = false
params.Mus_musculus = false
// Use with DESeq2 output only
// Vote Counting Default Method

process Join_Datasets { 
    conda "/media/run-projects/software/miniconda3/envs/metaVolcano" 
    publishDir "${params.outdir}/Data/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(organism_files)

    output:  
    path("${params.organism}.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript 
    library(MetaVolcanoR)

    datasets <- list()
    for (file in list.files("./", pattern = "_DE_NS.csv")){
        print(file)
        name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
        archivo <- read.delim(file,sep =",")
        colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
        archivo\$X <- NULL
        nam_best_grid_nosp <- name
        assign(nam_best_grid_nosp, archivo)
        #datasets <- append(datasets,get(nam_best_grid_nosp))
        datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)
    }
    saveRDS(datasets, file = "${params.organism}.rds")
    """   
}
process Join_Datasets_inter_specie { 
    conda "/media/run-projects/software/miniconda3/envs/metaVolcano" 
    publishDir "${params.outdir}/Data/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(organism_files)

    output:  
    path("${params.organism}.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript 
    library(MetaVolcanoR)

    datasets <- list()
    IDs_Hs <- read.delim("${params.eq_human}", sep= "\t")
    IDs_Mm <- read.delim("${params.eq_mouse}", sep= "\t")
    IDs_Mm\$gene_id <- sapply(strsplit(as.character(IDs_Mm\$gene_id),'.',fixed = TRUE), "[", 1)
    diopt <- read.delim("${params.diopt}", sep= "\t")[c("Search.Term","Ensmbl.ID...link.HPA.")]
    names(diopt) <- c("GeneID_Mm","GeneID_Hs")
    diopt <- merge(diopt,IDs_Hs,by.x="GeneID_Hs", by.y="gene_id")

    for (file in list.files("./", pattern = "_DE_NS.csv")){
        print(file)
        name <- gsub("_Control_vs_Case_DE_NS.csv","",file)
        archivo <- read.delim(file,sep =",")
        colnames(archivo)[colnames(archivo) =="Row.names"] <- "gene_id"
        archivo\$X <- NULL
        if (grepl("ENSMUSG",archivo\$gene_id[1])){
            archivo\$gene_name <- NULL
            archivo <- merge(archivo,diopt,by.x="gene_id",by.y="GeneID_Mm")
            archivo\$gene_id <- NULL
            names(archivo)[names(archivo)=="GeneID_Hs"] <- "gene_id"
            archivo <- archivo[order(archivo\$padj),]
            archivo <- archivo[!duplicated(archivo\$gene_id),]            
        }
        nam_best_grid_nosp <- name
        assign(nam_best_grid_nosp, archivo)
        datasets[[nam_best_grid_nosp]] <- get(nam_best_grid_nosp)
    }
    saveRDS(datasets, file = "${params.organism}.rds")
    """   
}
process Vote_counting { 
    conda "/media/run-projects/software/miniconda3/envs/metaVolcano" 
    publishDir "${params.outdir}", mode:'copy'
    maxForks 5   

    input:
    path(rds_file)

    output:  
    path("*")
    path("*_vote-count.rds")

    script:
    """
    #!/usr/bin/env Rscript 
    library(MetaVolcanoR)
    datasets <- readRDS("${rds_file}")
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


    important_genes <- meta_degs_vote@metaresult[meta_degs_vote@metaresult\$degvcount != "1.Unperturbed",]
    IDs_Hs <- read.delim("${params.eq_human}", sep= "\t")
    IDs_Mm <- read.delim("${params.eq_mouse}", sep= "\t")
    IDs_Mm\$gene_id <- sapply(strsplit(as.character(IDs_Mm\$gene_id),'.',fixed = TRUE), "[", 1)
    if (grepl("ENSMUSG",important_genes\$gene_id[1])){
            important_genes <- merge(IDs_Mm,important_genes,by="gene_id")
    } else {
            important_genes <- merge(IDs_Hs,important_genes,by="gene_id")
    }
    write.table(important_genes, file="${params.organism}_important_genes_votecounting_metaanalisis.tsv", row.names=F, sep="\t")
    saveRDS(important_genes, file = "${params.organism}_vote-count.rds")
    """   
}
process Random_Effect_Model { 
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
    library(MetaVolcanoR)
    datasets <- readRDS("${rds_file}")
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
    """   
}
process Combining_approach { 
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
    library(MetaVolcanoR)
    datasets <- readRDS("${rds_file}")
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
    """   
}
process Enrichr { 
    conda "/media/run-projects/software/miniconda3/envs/Trabajo" 
    publishDir "${params.outdir}", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(rds_file)

    output:  
    path("enrichment.pdf")

    script:
    """
    #!/usr/bin/env Rscript 
    library(enrichR)
    important_genes <- readRDS("${rds_file}")
    genes_alterados <- unique(important_genes\$gene_name)
    setEnrichrSite("Enrichr") 
    websiteLive <- TRUE
    dbs <- listEnrichrDbs()
    dbs <- c("$params.enrich_gene_set")
    enriched <- enrichr(genes_alterados, dbs)   
    enriched <- enriched[[1]][enriched[[1]]\$Adjusted.P.value < 0.05,]
    pdf("enrichment.pdf",width = 16, height=9)
    plotEnrich(enriched, showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value", title = paste(dbs[1], " enrichment\n", sep = "" ))
    dev.off()
    """   
}
workflow { 
    if (params.Inter_species != false){    
        Human = Channel.fromPath("${params.Dir}Homo_sapiens/*_DE_NS.csv")
        Mouse = Channel.fromPath("${params.Dir}Mus_musculus/*_DE_NS.csv")
        Organisms_ch = Human.concat(Mouse)
        params.organism = "All"
        params.eq_mouse = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.vM27.txt"
        params.eq_human = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.v38.txt"
        params.diopt = "/media/run-projects/Adolfo/MirScience/RNAseq/diopt8.5_results_Mm27_to_Hs.high.tab"
        params.enrich_gene_set= "GO_Biological_Process_2021"
    }
    else if (params.Mus_musculus != false){    
        Mouse = Channel.fromPath("${params.Dir}Mus_musculus/*_DE_NS.csv")
        Organisms_ch = Mouse
        params.organism = "Mus_musculus"
        params.eq_mouse = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.vM27.txt"
        params.eq_human = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.v38.txt"
        params.enrich_gene_set= "WikiPathways_2019_Mouse"
    } else {
        Human = Channel.fromPath("${params.Dir}Homo_sapiens/*_DE_NS.csv")
        Organisms_ch = Human
        params.organism = "Homo_sapiens"
        params.eq_mouse = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.vM27.txt"
        params.eq_human = "/media/run-projects/Adolfo/MirScience/RNAseq/IDs_equivalencias_gencode.v38.txt"
        params.enrich_gene_set= "GO_Biological_Process_2021"
    }       
    params.outdir = "Resultados/Meta-analysis/${params.organism}/${params.Term}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()
    if (params.Inter_species != false){ 
    Join_Datasets_inter_specie(Organisms_ch.collect())
    Vote_counting(Join_Datasets_inter_specie.out)
    } else {
    Join_Datasets(Organisms_ch.collect())
    Vote_counting(Join_Datasets.out)        
    }
    Enrichr(Vote_counting.out[1])
    //Dowload_Datasets(Runs)
}