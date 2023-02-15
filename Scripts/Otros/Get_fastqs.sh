#!/bin/bash
for R1 in $(cat SRR_Acc_List*.txt)
do
    echo $R1 
    fastq-dump --split-files --gzip $R1
done