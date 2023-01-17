#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.Organisms = ""
params.Disease_term = ""
params.All_organisms = false

process Dowload_Bioproject_info { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 5
    tag "Download datasets $accession"

    input:
    val(accession)

    output:  
    tuple val(accession),  val(GPL),  val(Bioproject), path("${accession}.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env python3
    import os
    from re import I
    from socket import if_indextoname
    import sys
    from os import walk
    from Bio import Entrez 
    import pandas as pd
    from bs4 import BeautifulSoup
    Entrez.email = "adolfo.rojas@ug.uchile.cl"
    Disease = "$params.Disease_term"
    db = 'bioproject'
    Query = Disease + '(' + Disease + '[Description] OR ' + Disease + '[Title]) AND ("transcriptome gene expression"[Filter] AND ("bioproject sra"[Filter] OR "bioproject gds"[Filter]))'
    res = Entrez.read(Entrez.esearch(db=db, term=Query,rettype="summary", retmode= "text", RetMax = 200000, usehistory = "y"))
    df = pd.DataFrame()        
    file = Entrez.read(Entrez.esummary(db=db, query_key = res["QueryKey"] , RetMax = 200000, retmode= "text", WebEnv = res["WebEnv"]))
    df2 = pd.DataFrame(file["DocumentSummarySet"]["DocumentSummary"])
    df = df.append(df2)
    print(df)
    df.to_excel("Bioprojects_"+Disease+".xlsx", index = False)
    df.to_csv("Bioprojects_"+Disease+".tsv",sep = "\t", index = False)
    """   
}
workflow { 
    params.outdir = "Resultados/${params.Disease_term}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()

    Sample_table = Channel.fromPath(params.GEO_TABLE)
    if (params.normalize == false){    
        print "\nCheck Information:\n $params.Project,$params.Organism\n"
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