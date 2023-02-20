#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.Organisms = ""
params.Disease_term = ""
params.All_organisms = false
// Mus musculus and Homo sapiens Default organisms 

process Dowload_Bioproject_info { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 5
    
    input:
    val(term)

    output:  
    tuple val(term), path("Bioprojects_full.tsv")

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
    Disease = "${term}"
    db = 'bioproject'
    Query = Disease + '(' + Disease + '[Description] OR ' + Disease + '[Title]) AND ("transcriptome gene expression"[Filter] AND ("bioproject sra"[Filter] OR "bioproject gds"[Filter]))'
    res = Entrez.read(Entrez.esearch(db=db, term=Query,rettype="summary", retmode= "text", RetMax = 200000, usehistory = "y"))
    df = pd.DataFrame()        
    file = Entrez.read(Entrez.esummary(db=db, query_key = res["QueryKey"] , RetMax = 200000, retmode= "text", WebEnv = res["WebEnv"]))
    df2 = pd.DataFrame(file["DocumentSummarySet"]["DocumentSummary"])
    df = df.append(df2)
    print(df)
    #df.to_excel("Bioprojects_full.xlsx", index = False)
    df.to_csv("Bioprojects_full.tsv",sep = "\t", index = False)
    """   
}
process Dowload_SRA_info { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 10
    tag "Download datasets $accession"

    input:
    tuple val(term), path(Bioproject_TSV)

    output:  
    tuple val(term), path("*_Biosample_SRA.tsv")

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
    Disease = "${term}"
    df = pd.read_csv("${Bioproject_TSV}", sep="\t")

    db = 'bioproject'
    df5 = pd.DataFrame()
    not_in_SRA = []
    df4 = pd.DataFrame()
    for SRA_acc in df.loc[(df["TaxId"] == 10090) | (df["TaxId"] == 9606)]["Project_Acc"].unique():
        try:
            Query3 = SRA_acc + '[BioProject] AND ("biomol rna"[Properties] AND "filetype fastq"[Properties])'
            res3 = Entrez.read(Entrez.esearch(db="sra", term=Query3, retmode= "xml", RetMax = 200000, usehistory = "y"))                 
            file2 = Entrez.read(Entrez.esummary(db="sra", query_key = res3["QueryKey"], RetMax = 200000,WebEnv = res3["WebEnv"]))    
            df3_in_loop = pd.DataFrame(file2)
            df3_in_loop["BioProjectID"] = SRA_acc
            if len(df3_in_loop["Id"].unique()) >= 6:         
                df4 = df4.append(df3_in_loop)
                biosam = Entrez.read(Entrez.elink(dbfrom="sra",db="biosample", linkname="sra_biosample", WebEnv = res3["WebEnv"],query_key = res3["QueryKey"]))            
                l = biosam[0]["LinkSetDb"][0]["Link"]
                biosam_ids = [d['Id'] for d in l]
                file = Entrez.read(Entrez.esummary(db="biosample",id = ",".join(biosam_ids),RetMax = 2000, retmode= "xml"))               
                df4_in_loop = pd.DataFrame(file["DocumentSummarySet"]["DocumentSummary"]).drop(columns={"Package","SortKey","PublicationDate","Date","ModificationDate","Infraspecies","Identifiers"})               
                df5 = df5.append(df4_in_loop)
                print(df5)
        except RuntimeError:
            not_in_SRA.append(SRA_acc)
            pass

    df4["Run_acc"] = df4["Runs"].str.split('Run acc="',expand=True)[1].str.split('"', expand=True)[0] 
    df4["is_public"] = df4["Runs"].str.split('is_public="',expand=True)[1].str.split('"', expand=True)[0] 

    df4["Library_name"] = df4["ExpXml"].str.split("LIBRARY_NAME>",expand=True)[1].str.split("</", expand=True)[0]
    df4["Library_Strategy"] = df4["ExpXml"].str.split("IBRARY_STRATEGY>",expand=True)[1].str.split("</", expand=True)[0]   
    df4["Library_Source"] = df4["ExpXml"].str.split("LIBRARY_SOURCE>",expand=True)[1].str.split("</", expand=True)[0]   
    df4["Library_Selection"] = df4["ExpXml"].str.split("LIBRARY_SELECTION>",expand=True)[1].str.split("</", expand=True)[0]   
    df4["Library_Layout"] = df4["ExpXml"].str.split("LIBRARY_LAYOUT>",expand=True)[1].str.split("</", expand=True)[0].str.replace("<","").str.replace("/>","")   
    df4["Library_Construction_Protocol"] = df4["ExpXml"].str.split("LIBRARY_CONSTRUCTION_PROTOCOL>",expand=True)[1].str.split("</", expand=True)[0]  
    df4["Title"] = df4["ExpXml"].str.split("Title>",expand=True)[1].str.split("</", expand=True)[0]  
    df4["BiosampleID"] = df4["ExpXml"].str.split("Biosample>",expand=True)[1].str.split("</", expand=True)[0]   
    df4["Organism_name"] = df4["ExpXml"].str.split('ScientificName="',expand=True)[1].str.split('"', expand=True)[0]  
    df4["Platform"] = df4["ExpXml"].str.split('Platform instrument_model="',expand=True)[1].str.split('"', expand=True)[0]  
    df4["Study_name"] = df4["ExpXml"].str.split(r'><Study \\S+ name="',expand=True)[1].str.split('"', expand=True)[0] 
    df4 = df4.drop(columns={"ExpXml","Runs","Item","ExtLinks","Id","Library_name","UpdateDate","CreateDate"})

    df5 = df5.reset_index()
    contador = 0
    df5_attr = dict()
    df5_attr["Accession"] = []
    for i in range(0,len(df5["SampleData"])):
        S = BeautifulSoup(df5["SampleData"][i], 'lxml')
        attr_prsentes = []
        attr_prsentes.append("Accession")
        df5_attr["Accession"].append(df5["Accession"][i])    
        for ii in S.find_all("attribute"):
            attr_prsentes.append(ii.get('attribute_name'))        
            if contador == 0:
                df5_attr[ii.get('attribute_name')] = [ii.text]
            else:
                if ii.get('attribute_name') not in df5_attr.keys():
                    df5_attr[ii.get('attribute_name')] = [""]* contador
                df5_attr[ii.get('attribute_name')].append(ii.text)
        for iii  in [x for x in df5_attr.keys() if x not in attr_prsentes]:
            df5_attr[iii].append("")
        contador+=1

    df5_attr_df = pd.DataFrame(df5_attr)
    column_order = dict()
    column_order["column"]= []
    column_order["non_NA_values"]= []
    for i in df5_attr_df.columns:
        column_order["non_NA_values"].append(len(df5_attr_df[df5_attr_df[i].str.contains(r'\\w')].index))
        column_order["column"].append(i)
    column_order_df = pd.DataFrame(column_order).sort_values(by=['non_NA_values'], ascending=False)
    df5_attr_df = df5_attr_df[column_order_df["column"]]

    samplesNumber = df4.groupby("BioProjectID")["BiosampleID"].count().to_frame().rename(columns={"BiosampleID":"samplesNumber"})
    df4 = pd.merge(df4,samplesNumber,left_on="BioProjectID",right_index=True)
    df5 = df5.drop(columns={"Taxonomy","SourceSample","SampleData"})
    df5_backup = df5.copy()
    df5 = df5.merge(df5_attr_df, on="Accession")
    df5.to_excel(Disease.replace(" ","_") +"_Biosample.xlsx", index = False)
    df5.to_csv(Disease.replace(" ","_")+"_Biosample.tsv",sep = "\t", index = False)

    df8 = df4.merge(df5,left_on="BiosampleID",right_on="Accession")
    df8.to_excel(Disease.replace(" ","_")+"_Biosample_SRA.xlsx", index = False)
    df8.to_csv(Disease.replace(" ","_")+"_Biosample_SRA.tsv",sep = "\t", index = False)
    df8 = df8[["BioProjectID","Run_acc","Library_Layout","Platform","Organism"]]
    df8["Condition"] = "ND"
    df8["Organism"] = df8["Organism"].str.replace(" ","_")
    df8["Library_Layout"] = df8["Library_Layout"].str.replace(" ","")
    df8["Platform"] = df8["Platform"].str.replace(" ","_").str.lower()
    df8.to_csv(Disease.replace(" ","_")+"_samples_SRA.csv",sep = ",", index = False)
    """   
}
process Dowload_SRA_info_v2 { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/SRA/${Bioproject}", mode:'copy', pattern: '*.tsv' 
    maxForks 7   
    errorStrategy 'retry'
    maxRetries 5
    tag "Download datasets $Bioproject"

    input:
    val(Bioproject)

    output:  
    tuple val(Bioproject), path("*.tsv"), optional: true, emit: SRA
    tuple val(Bioproject), path("*.txt"), optional: true, emit: not_SRA

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
    SRA_acc = "${Bioproject}"

    df5 = pd.DataFrame()
    not_in_SRA = []
    df4 = pd.DataFrame()
    try:
        Query3 = SRA_acc + '[BioProject] AND ("biomol rna"[Properties] AND "filetype fastq"[Properties])'
        res3 = Entrez.read(Entrez.esearch(db="sra", term=Query3, retmode= "xml", RetMax = 200000, usehistory = "y"))                 
        file2 = Entrez.read(Entrez.esummary(db="sra", query_key = res3["QueryKey"], RetMax = 200000,WebEnv = res3["WebEnv"]))    
        df3_in_loop = pd.DataFrame(file2)
        df3_in_loop["BioProjectID"] = SRA_acc
        if len(df3_in_loop["Id"].unique()) >= 6:         
            df4 = df4.append(df3_in_loop)
            biosam = Entrez.read(Entrez.elink(dbfrom="sra",db="biosample", linkname="sra_biosample", WebEnv = res3["WebEnv"],query_key = res3["QueryKey"]))            
            l = biosam[0]["LinkSetDb"][0]["Link"]
            biosam_ids = [d['Id'] for d in l]
            file = Entrez.read(Entrez.esummary(db="biosample",id = ",".join(biosam_ids),RetMax = 2000, retmode= "xml"))               
            df4_in_loop = pd.DataFrame(file["DocumentSummarySet"]["DocumentSummary"]).drop(columns={"Package","SortKey","PublicationDate","Date","ModificationDate","Infraspecies","Identifiers"})               
            df5 = df5.append(df4_in_loop)
            print(df5)
            df4["Run_acc"] = df4["Runs"].str.split('Run acc="',expand=True)[1].str.split('"', expand=True)[0] 
            df4["is_public"] = df4["Runs"].str.split('is_public="',expand=True)[1].str.split('"', expand=True)[0] 
            df4["Library_Strategy"] = df4["ExpXml"].str.split("IBRARY_STRATEGY>",expand=True)[1].str.split("</", expand=True)[0]   
            df4["Library_Source"] = df4["ExpXml"].str.split("LIBRARY_SOURCE>",expand=True)[1].str.split("</", expand=True)[0]   
            df4["Library_Selection"] = df4["ExpXml"].str.split("LIBRARY_SELECTION>",expand=True)[1].str.split("</", expand=True)[0]   
            df4["Library_Layout"] = df4["ExpXml"].str.split("LIBRARY_LAYOUT>",expand=True)[1].str.split("</", expand=True)[0].str.replace("<","").str.replace("/>","")  
            try: 
                df4["Library_Construction_Protocol"] = df4["ExpXml"].str.split("LIBRARY_CONSTRUCTION_PROTOCOL>",expand=True)[1].str.split("</", expand=True)[0] 
            except KeyError:
                print("error") 
            df4["Title"] = df4["ExpXml"].str.split("Title>",expand=True)[1].str.split("</", expand=True)[0]  
            df4["BiosampleID"] = df4["ExpXml"].str.split("Biosample>",expand=True)[1].str.split("</", expand=True)[0]   
            df4["Organism_name"] = df4["ExpXml"].str.split('ScientificName="',expand=True)[1].str.split('"', expand=True)[0]  
            df4["Platform"] = df4["ExpXml"].str.split('Platform instrument_model="',expand=True)[1].str.split('"', expand=True)[0]  
            df4["Study_name"] = df4["ExpXml"].str.split(r'><Study \\S+ name="',expand=True)[1].str.split('"', expand=True)[0] 
            df4 = df4.drop(columns={"ExpXml","Runs","Item","ExtLinks","Id","UpdateDate","CreateDate"})
            df5 = df5.reset_index()
            contador = 0
            df5_attr = dict()
            df5_attr["Accession"] = []
            for i in range(0,len(df5["SampleData"])):
                S = BeautifulSoup(df5["SampleData"][i], 'lxml')
                df5_attr["Accession"].append(df5["Accession"][i])
                df5_attr2 = dict() 
                df5_attr2["Accession"] = []
                df5_attr2["Accession"].append(df5["Accession"][i])
                for ii in S.find_all("attribute"):   
                    if ii.get('attribute_name') in df5_attr2.keys():
                        df5_attr2[ii.get('attribute_name')+"_2"] = [] 
                        if contador == 0:
                            df5_attr[ii.get('attribute_name')+"_2"] = [ii.text]
                            df5_attr2[ii.get('attribute_name')+"_2"] = [ii.text]
                        else:
                            df5_attr[ii.get('attribute_name')+"_2"].append(ii.text)   
                            df5_attr2[ii.get('attribute_name')+"_2"].append(ii.text)   
                    else:   
                        df5_attr2[ii.get('attribute_name')] = []       
                        if contador == 0:
                            df5_attr[ii.get('attribute_name')] = [ii.text]
                            df5_attr2[ii.get('attribute_name')] = [ii.text]
                        else:
                            df5_attr[ii.get('attribute_name')].append(ii.text)
                            df5_attr2[ii.get('attribute_name')].append(ii.text)
                contador+=1
            df5_attr_df = pd.DataFrame(df5_attr)
            column_order = dict()
            column_order["column"]= []
            column_order["non_NA_values"]= []
            for i in df5_attr_df.columns:
                column_order["non_NA_values"].append(len(df5_attr_df[df5_attr_df[i].str.contains(r'\\w')].index))
                column_order["column"].append(i)
            column_order_df = pd.DataFrame(column_order).sort_values(by=['non_NA_values'], ascending=False)
            df5_attr_df = df5_attr_df[column_order_df["column"]]
            samplesNumber = df4.groupby("BioProjectID")["BiosampleID"].count().to_frame().rename(columns={"BiosampleID":"samplesNumber"})
            df4 = pd.merge(df4,samplesNumber,left_on="BioProjectID",right_index=True)
            df5 = df5.drop(columns={"Taxonomy","SourceSample","SampleData"})
            df5_backup = df5.copy()
            df5 = df5.merge(df5_attr_df, on="Accession")
            df8 = df4.merge(df5,left_on="BiosampleID",right_on="Accession")
            df8.to_csv("Biosample_SRA.tsv",sep = "\t", index = False)
            df8 = df8[["BioProjectID","Run_acc","Library_Layout","Platform","Organism"]]
            df8["Condition"] = "ND"
            df8["Organism"] = df8["Organism"].str.replace(" ","_")
            df8["Library_Layout"] = df8["Library_Layout"].str.replace(" ","")
            df8["Platform"] = df8["Platform"].str.replace(" ","_").str.lower()
            df8.to_csv("Samples_SRA.tsv",sep = "\t", index = False)
    except RuntimeError:
        not_in_SRA = pd.DataFrame({"Project": ["${Bioproject}"]})
        not_in_SRA.to_csv("Not_SRA.txt",sep = "\t", index = False)

    """   
}
process Dowload_GEO_info { 
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
    not_in_GEO = []
    df3 = pd.DataFrame()
    muestras_GEO_affi = pd.DataFrame()
    for geo_acc in df.loc[df["Project_Acc"].isin(not_in_SRA)]["Project_Acc"].unique():
        try:
            Query2 = geo_acc + '[All Fields] AND ("Expression profiling by array"[Filter] OR "Expression profiling by high throughput sequencing"[Filter] OR "Expression profiling by RT-PCR"[Filter] OR "Non-coding RNA profiling by high throughput sequencing"[Filter] OR "Non-coding RNA profiling by genome tiling array"[Filter] OR "Expression profiling by genome tiling array"[Filter] OR "Non-coding RNA profiling by array"[Filter] OR "Expression profiling by SNP array"[Filter] OR "Expression profiling by SAGE"[Filter] OR "Expression profiling by MPSS"[Filter])'
            res2 = Entrez.read(Entrez.esearch(db="gds", term=Query2,rettype="summary", retmode= "xml", RetMax = 200000, usehistory = "y"))               
            file = Entrez.read(Entrez.esummary(db="gds", query_key = res2["QueryKey"], RetMax = 2000, retmode= "xml", WebEnv = res2["WebEnv"]))              
            df2_in_loop = pd.DataFrame(file)
            df2_in_loop["BioProjectID"] = geo_acc
            if len(df2_in_loop.loc[df2_in_loop["n_samples"]>= 6]) > 0:
                df3 = df3.append(df2_in_loop)
                print(df3)
                muestras_GEO_affi_loop = pd.DataFrame(file[0]["Samples"])
                muestras_GEO_affi_loop["Accession2"] = df2_in_loop["Accession"][0]
                muestras_GEO_affi_loop["GPL"] = df2_in_loop["GPL"][0]
                muestras_GEO_affi_loop["BioProjectID"] = df2_in_loop["BioProjectID"][0]
                muestras_GEO_affi = muestras_GEO_affi.append(muestras_GEO_affi_loop)
        except RuntimeError:
            not_in_GEO.append(geo_acc)
            pass

    muestras_GEO_affi["Condition"] = "ND"
    df3 = df3.loc[df3["entryType"]== "GSE"]                
    df3.to_excel("MirScience/Datasets/BioProject/"+Disease.replace(" ","_")+"_GEO_projects.xlsx", index = False)
    df3.to_csv("MirScience/Datasets/BioProject/"+Disease.replace(" ","_")+"_GEO_projects.tsv",sep = "\t", index = False)
    df3 = df3[["BioProjectID","taxon"]]
    df3["taxon"] = df3["taxon"].str.replace(" ","_")
    muestras_GEO_affi = pd.merge(muestras_GEO_affi,df3,on="BioProjectID")
    muestras_GEO_affi = muestras_GEO_affi[["Accession","Title","Accession2","GPL","BioProjectID","taxon","Condition"]]
    muestras_GEO_affi.to_csv("MirScience/Datasets/BioProject/"+Disease.replace(" ","_")+"_GEO_samples.csv",sep = ",", index = False)
    df_not_GEO = pd.DataFrame(not_in_GEO)
    df_not_GEO.to_csv("MirScience/Datasets/BioProject/"+Disease.replace(" ","_")+"_not_found_projects.txt")
    """   
}
process Dowload_GEO_info_v2 { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/GEO/${Bioproject}", mode:'copy', pattern: '*.tsv' 
    maxForks 3   
    errorStrategy 'retry'
    maxRetries 5
    tag "Download datasets $Bioproject"

    input:
    val(Bioproject)

    output:  
    tuple val(Bioproject), path("*.tsv"), optional: true, emit: GEO
    tuple val(Bioproject), path("*.txt"), optional: true, emit: not_GEO

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
    geo_acc = "${Bioproject}"
    not_in_GEO = []
    df3 = pd.DataFrame()
    muestras_GEO_affi = pd.DataFrame()
    try:
        Query2 = geo_acc + '[All Fields] AND ("Expression profiling by array"[Filter] OR "Expression profiling by high throughput sequencing"[Filter] OR "Expression profiling by RT-PCR"[Filter] OR "Non-coding RNA profiling by high throughput sequencing"[Filter] OR "Non-coding RNA profiling by genome tiling array"[Filter] OR "Expression profiling by genome tiling array"[Filter] OR "Non-coding RNA profiling by array"[Filter] OR "Expression profiling by SNP array"[Filter] OR "Expression profiling by SAGE"[Filter] OR "Expression profiling by MPSS"[Filter])'
        res2 = Entrez.read(Entrez.esearch(db="gds", term=Query2,rettype="summary", retmode= "xml", RetMax = 200000, usehistory = "y"))               
        file = Entrez.read(Entrez.esummary(db="gds", query_key = res2["QueryKey"], RetMax = 2000, retmode= "xml", WebEnv = res2["WebEnv"]))              
        df2_in_loop = pd.DataFrame(file)
        df2_in_loop["BioProjectID"] = geo_acc
        if len(df2_in_loop.loc[df2_in_loop["n_samples"]>= 6]) > 0:
            df3 = df3.append(df2_in_loop)
            print(df3)
            muestras_GEO_affi_loop = pd.DataFrame(file[0]["Samples"])
            muestras_GEO_affi_loop["Accession2"] = df2_in_loop["Accession"][0]
            muestras_GEO_affi_loop["GPL"] = df2_in_loop["GPL"][0]
            muestras_GEO_affi_loop["BioProjectID"] = df2_in_loop["BioProjectID"][0]
            muestras_GEO_affi = muestras_GEO_affi.append(muestras_GEO_affi_loop)
            muestras_GEO_affi["Condition"] = "ND"
            df3 = df3.loc[df3["entryType"]== "GSE"]                
            df3.to_csv("GEO_projects.tsv",sep = "\t", index = False)
            df3 = df3[["BioProjectID","taxon"]]
            df3["taxon"] = df3["taxon"].str.replace(" ","_")
            muestras_GEO_affi = pd.merge(muestras_GEO_affi,df3,on="BioProjectID")
            muestras_GEO_affi = muestras_GEO_affi[["Accession","Title","Accession2","GPL","BioProjectID","taxon","Condition"]]
            muestras_GEO_affi.to_csv("GEO_samples.tsv",sep = "\t", index = False)
    except RuntimeError:
        not_in_GEO = pd.DataFrame({"Project": ["${Bioproject}"]})
        not_in_GEO.to_csv("Not_GEO.txt",sep = "\t", index = False)
    """   
}
process Clean_output_dir_SRA { 
    maxForks 10   

    input:
    val(Bioproject)

    output:
    stdout

    script:
    """
    if [ `ls ../../../${params.outdir}/SRA/${Bioproject}/Biosample_SRA.tsv > /dev/null ; echo \$?` -eq 0 ]; then ls -lh ../../../${params.outdir}/SRA/${Bioproject}/Biosample_SRA.tsv ; else rm -r ../../../${params.outdir}/SRA/${Bioproject}/; fi
    """   
}
process Clean_output_dir_GEO { 
    maxForks 10   

    input:
    val(Bioproject)

    output:
    stdout

    script:
    """
    if [ `ls ../../../${params.outdir}/GEO/${Bioproject}/GEO_samples.tsv > /dev/null ; echo \$?` -eq 0 ]; then ls -lh ../../../${params.outdir}/GEO/${Bioproject}/GEO_samples.tsv ; else rm -r ../../../${params.outdir}/GEO/${Bioproject}/; fi
    """   
}
process Refine_Bioproject { 
    conda "/media/run-projects/software/miniconda3/envs/Personal" 
    publishDir "${params.outdir}/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 1
    beforeScript "echo ${Projectos_validos} | sed 's/, /\\n/g' | sed 's/[\\[]//g' | sed 's/[]]//g' > projectos_validos.txt"
    
    input:
    tuple val(term), path("Bioprojects_full_tsv")
    val(Projectos_validos)


    output:  
    path("Bioprojects_Reduced.tsv")

    script:
    """
    #!/usr/bin/env python3
    import os
    import sys
    from os import walk
    import pandas as pd
    Datos = pd.read_csv("$Bioprojects_full_tsv",sep="\t")
    Lita_projectos_validos = pd.read_csv("projectos_validos.txt",sep="\t",header=None)
    Datos = Datos.loc[Datos["Project_Acc"].isin(Lita_projectos_validos[0])]
    Datos.to_csv("Bioprojects_Reduced.tsv",sep = "\t", index = False)
    """   
}
workflow { 
    Directorio = params.Disease_term.replaceAll(/ /, '_')
    params.outdir = "Resultados/Datasets/${Directorio}"
    out_dir = file(params.outdir)
    out_dir.mkdir()    

    print 'If something is wrong press "ctrl + C" and change parameters'
    print "Workflow project outdir " + params.outdir
    Dowload_Bioproject_info(params.Disease_term)
    ST_human = Dowload_Bioproject_info.out.map{it[1]}.splitCsv(header: true,sep:"\t").filter(row -> row.TaxId == "9606")
    ST_mouse = Dowload_Bioproject_info.out.map{it[1]}.splitCsv(header: true,sep:"\t").filter(row -> row.TaxId == "10090")
    ST = ST_mouse.concat(ST_human)
    ST = ST.map { row-> row.Project_Acc}
    Dowload_SRA_info_v2(ST)
    Dowload_GEO_info_v2(Dowload_SRA_info_v2.out.not_SRA.map{it[0]})
    Available_projects = Dowload_GEO_info_v2.out.GEO.map{it[0]}.concat(Dowload_SRA_info_v2.out.SRA.map{it[0]}).collect()
    Refine_Bioproject(Dowload_Bioproject_info.out,Available_projects)
    Clean_output_dir_SRA(Dowload_SRA_info_v2.out.not_SRA.map{it[0]})
    Clean_output_dir_GEO(Dowload_GEO_info_v2.out.not_GEO.map{it[0]})
}