cd /media/storage/Adolfo/MirScience/RNA-seq/prueba_short_reads/miRNAs/
fastqc *
multiqc *.zip -n reporte
scp -P 1313 reporte.html adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/reporte_miRNAseq.html

for run in $(cat ../Acc_List.txt)
do 
fastp --thread 60 --in1 raw/${run}_1.fastq.gz --length_required 18 --length_limit 30 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --out1 trimmedV2/${run}_1_trimmed_R1.fastq.gz
done

fastqc *
multiqc *.zip -n reporte_miRNAseq2
scp -P 1313 reporte_miRNAseq2.html adolfo@200.89.65.156:/media/run-projects/Adolfo/MirScience/

bowtie-build /media/storage/datasets/genomes/mouse/Mm_GRCm39/GRCm39.genome.fa /media/storage/datasets/indexes/bowtie/mouse/Mm_GRCm39/GRCm39.genome.fa
extract_miRNAs.pl hairpin.fa.gz mmu > hairpin_ref.fa
extract_miRNAs.pl hairpin.fa mmu > hairpin_ref.fa

mapper.pl SRR18454579_1_trimmed_R1.fastq -e -h -i -j -l 18 -m -p genome/GRCm39.genome.fa -s reads_collapsed.fa -t reads_vs_refdb.arf -v -o 60

miRDeep2.pl reads_collapsed.fa genome/GRCm39.genome.fa reads_vs_refdb.arf mature_ref.fa none hairpin_ref.fa -t mmu >report.log


#### quantifier

gunzip -c ../../trimmedV2/SRR18454578_1_trimmed_R1.fastq.gz > SRR18454578_1_trimmed_R1.fastq

fastq2fasta.pl SRR18454578_1_trimmed_R1.fastq > SRR18454578.fa

collapse_reads_md.pl SRR18454578.fa mmu > SRR18454578.collapsed.fa
quantifier.pl -p ../hairpin_ref.fa -m ../mature_ref.fa -r SRR18454578.collapsed.fa -t mmu -d -j




import pandas as pd 
Datos = pd.read_csv("diopt8.5_results_Mm27_to_Hs.csv")[['Search Term', 'Mouse GeneID', 'Mouse Symbol','Human GeneID', 'Human Species Gene ID', 'Human Symbol','Ensmbl ID  (link HPA)', 'DIOPT Score', 'Weighted Score', 'Rank']]
Datos = Datos.loc[Datos["Rank"]!="low"]
Datos = Datos[Datos['Rank'].notna()]
Datos = Datos[Datos['Ensmbl ID  (link HPA)'].notna()]
Datos = Datos.loc[Datos["Rank"]!="moderate"]
Datos.to_csv("diopt8.5_results_Mm27_to_Hs.high.tab",sep="\t", index=False)


#sample_df = sample_df.groupby(by=sample_df.columns, axis=1).sum()