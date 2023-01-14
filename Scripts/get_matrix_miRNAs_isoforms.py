#!/usr/bin/env python3
import os
import pandas as pd

directorios = list()
for root, dirs, files in os.walk(".", topdown=False):
    for name in dirs:
        dir = os.path.join(root, name)
        directorios.append(dir)

final_df = pd.DataFrame()
count = 0
info  = pd.read_csv("IDs_TCGA_samples.tab", sep="\t", header=0)
for FILE in directorios:
        for root, dirs, files in os.walk(FILE, topdown=False):
                for name in files:
                        count+= 1                        
                        file = os.path.join(root, name)
                        archivo = pd.read_csv(file, sep="\t", header= 0)
                        archivo = archivo[archivo.miRNA_region.str.contains('mature,MIMAT', na=False)]
                        archivo["miRNA_region"] = archivo["miRNA_region"].str.split(',', n=-1, expand = True)[1]                         
                        data = {}
                        for mirna_mature in archivo["miRNA_region"].unique():
                            data[mirna_mature] = archivo[archivo.miRNA_region == mirna_mature].read_count.sum()                        
                        file = file.split('/')[2]                                              
                        df = pd.DataFrame.from_dict(data, orient='index', columns=info[info.file_name == file].cases)
                        final_df = df.join(final_df, how = "outer").fillna(0)
                        print(final_df)
                        print(count)

final_df.to_csv("matrix_miRNAs_isoforms.tab", sep="\t",header=True, index=True)
