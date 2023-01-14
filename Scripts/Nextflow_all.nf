#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Mandatory Params
params.reads = "Raw_data/*_{1,2}.fastq.gz" // Raw reads a analizar
params.outdir = "Resultados" // Directorio donde se guardaran archivos de interes
params.ref_genome = "output_data/genome.fasta"  // Genoma de referencia de la especie
// Optional Params
params.ploidy = "1" // ploidia del organismo analizado
params.update_joint_vcf = true // existe un analisis previo de llamado de variantes
params.previus_vcf = "output_data/*.HC.raw.g.vcf.gz*" // escribir directorio de base de datos previa para actualizar o de los vcfs, si se desea crear una base de datos nueva
params.ref_seq_header = "NC_010168.1 Renibacterium salmoninarum ATCC 33209, complete sequence"
params.ref_seq_new_header ="ATCC_33209"
params.phame_conf_file = "renibacteriumsalmoninarum.ctl"
params.Trimmed = "${params.outdir}/Trimmed/*_trimmed_R{1,2}.fastq.gz"
// Channels definition
ref_genome_ch = channel.fromPath(params.ref_genome)
read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true) 
out_dir = file(params.outdir)
out_dir.mkdir()

process FASTQC {
    publishDir "${params.outdir}/FASTQC", mode:'copy', pattern: '*.html'      
    tag "FASTQC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:    
    path "*.html"
    path "*zip", emit: multiqc_input

    script:
    """
    fastqc $reads
    """
}
process MULTIQC {    
    publishDir params.outdir, mode:'copy'  

    input:
    file archivos
    val(name)

    output:
    path "*_multiqc_report.html"
    
    script:
    """    
    multiqc $archivos -n '$name'_multiqc_report.html
    """
}
process FASTP {    
    maxForks 5    
    publishDir "${params.outdir}/Trimmed", mode:'copy', pattern: '*.fastq.gz' 
    tag "FASTP on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads      
        
    script:
    """
    fastp --thread $task.cpus --in1 ${reads[0]} --in2 ${reads[1]} --length_required 40 --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_tail1=5 --trim_tail2=5 --trim_front1=15 --trim_front2=15 --out1 '$sample_id'_trimmed_R1.fastq.gz --out2 '$sample_id'_trimmed_R2.fastq.gz
    """
}
process REF_INDEX {            
    
    input:
    path x

    output:
    path("${x}.*"), emit: indexs
    path "*.dict", emit: dict

    script:
    """    
    bwa index -a is $x -p $x
    samtools faidx $x
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CreateSequenceDictionary -R $x
    java -jar /media/storage2/software/Picard/picard.jar ScatterIntervalsByNs REFERENCE= $x OUTPUT= '$x'.interval_list
    """
}
process BWA_MEM {    
    maxForks 5
    cpus 30
    tag "BWA_MEM on $sample_id"       
    
    input:
    tuple val(sample_id), path(reads)
    each path(index_files)
    each path(x)
    each path(y)

    output:
    path("*HC.raw.g.vcf.gz*"), emit: raw_g_vcf
        
    script:
    """
    bwa mem -v 1 -M -t $task.cpus $x ${reads[0]} ${reads[1]} > '$sample_id'_aligned.sam
    samtools view -S -b '$sample_id'_aligned.sam > '$sample_id'_aligned.bam
    java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'_aligned.bam OUTPUT= '$sample_id'_aligned_sorted.bam
    java -jar /media/storage2/software/Picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT I= '$sample_id'_aligned_sorted.bam O= '$sample_id'_aligned_dups_removed.bam REMOVE_DUPLICATES=true M=metrics
    java -jar /media/storage2/software/Picard/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I= '$sample_id'_aligned_dups_removed.bam O= '$sample_id'_RG.bam SO=coordinate RGLB=lib_1 RGPL=illumina RGPU=barcode_1 RGSM='$sample_id' CREATE_INDEX=true
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar HaplotypeCaller -R $x -I '$sample_id'_RG.bam -O '$sample_id'.HC.raw.g.vcf.gz -ERC GVCF --verbosity WARNING --native-pair-hmm-threads $task.cpus --sample-ploidy $params.ploidy
    rm *aligned* 
    """
}
process DBImport {    
    maxForks 5
    cpus 30
        
    input:
    path(raw_g_vcfs)
    path(index_files)
    path(x)
    path(y)

    output:
    tuple val("${raw_g_vcfs.size()/2}"), path("*.GenotypeGVCFs.raw.vcf.gz*"), emit: GenotypeGVCFs_raw  
        
    script:
    """
    ls *.HC.raw.g.vcf.gz > gvcf.list
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenomicsDBImport --genomicsdb-workspace-path GenomicsDB/ --intervals '$x'.interval_list --variant gvcf.list  
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenotypeGVCFs -R '$x' -V gendb://GenomicsDB/ --sample-ploidy $params.ploidy -O ${raw_g_vcfs.size()/2}.GenotypeGVCFs.raw.vcf.gz
    """
}
process DBImport_update_prev_vcf {    
    maxForks 5
    cpus 30
    // utiliza vcf previos y genera una base de datos nueva
        
    input:
    path(raw_g_vcfs)
    path(index_files)
    path(x)
    path(y)   
    path(previus_vcf) 

    output:    
    tuple val("${(raw_g_vcfs.size() + previus_vcf.size())/2}"), path("*.GenotypeGVCFs.raw.vcf.gz*"), emit: GenotypeGVCFs_raw  
        
    script:
    """
    ls *.HC.raw.g.vcf.gz > gvcf.list
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenomicsDBImport --genomicsdb-workspace-path GenomicsDB/ --intervals '$x'.interval_list --variant gvcf.list  
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenotypeGVCFs -R '$x' -V gendb://GenomicsDB/ --sample-ploidy $params.ploidy -O ${(raw_g_vcfs.size() + previus_vcf.size())/2}.GenotypeGVCFs.raw.vcf.gz
    """
}
process DBImport_update_db {    
    maxForks 5
    cpus 30
    // actualiza una base de datos preexistente
        
    input:
    path(raw_g_vcfs)
    path(index_files)
    path(x)
    path(y)   
    path(previus_db) 

    output:
    tuple val("${raw_g_vcfs.size()/2}_added"), path("*.GenotypeGVCFs.raw.vcf.gz*"), emit: GenotypeGVCFs_raw    
        
    script:
    """
    ls *.HC.raw.g.vcf.gz > gvcf.list
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenomicsDBImport --genomicsdb-update-workspace-path '$previus_db' --intervals '$x'.interval_list --variant gvcf.list  
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenotypeGVCFs -R '$x' -V gendb://'$previus_db' --sample-ploidy $params.ploidy -O ${raw_g_vcfs.size()/2}_added.GenotypeGVCFs.raw.vcf.gz
    """
}
process VCF_filters {    
    maxForks 5
    // actualiza una base de datos preexistente
        
    input:
    tuple val(name), path(joint_vcf)
    path(genoma)    
    path(index_files)
    path(dict)

    output:
    path "*.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz*", emit: vcf    
    path "GVCF_sample_list", emit: sample_list    
        
    script:
    """
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R '$genoma' -V ${joint_vcf[0]} -select-type SNP -O '$name'.GenotypeGVCFs.raw.SNPs.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantFiltration -R '$genoma' -V '$name'.GenotypeGVCFs.raw.SNPs.vcf.gz -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O '$name'.GenotypeGVCFs.hard_filtered.SNPs.vcf.gz    
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R '$genoma' -V ${joint_vcf[0]} -select-type INDEL -O '$name'.GenotypeGVCFs.raw.INDELs.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantFiltration -R '$genoma' -V '$name'.GenotypeGVCFs.raw.INDELs.vcf.gz -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' -O '$name'.GenotypeGVCFs.hard_filtered.INDELs.vcf.gz
    vcf-concat '$name'.GenotypeGVCFs.hard_filtered.SNPs.vcf.gz '$name'.GenotypeGVCFs.hard_filtered.INDELs.vcf.gz > '$name'.GenotypeGVCFs.hard_filtered.all.vcf
    vcftools --vcf '$name'.GenotypeGVCFs.hard_filtered.all.vcf --recode --stdout --remove-filtered-all | vcf-sort  > '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf
    bgzip -c '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf > '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz
    tabix -p vcf '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz
    vcf-query -l '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz > GVCF_sample_list
    """
}
process VCF_filters_GATK_only {    
    maxForks 5
    // actualiza una base de datos preexistente
        
    input:
    tuple val(name), path(joint_vcf)
    path(genoma)    
    path(index_files)
    path(dict)

    output:
    path "*.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz*", emit: vcf    
    path "GVCF_sample_list", emit: sample_list    
        
    script:
    """
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R '$genoma' -V ${joint_vcf[0]} -select-type SNP -O '$name'.GenotypeGVCFs.raw.SNPs.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantFiltration -R '$genoma' -V '$name'.GenotypeGVCFs.raw.SNPs.vcf.gz -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O '$name'.GenotypeGVCFs.hard_filtered.SNPs.vcf.gz    
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R '$genoma' -V ${joint_vcf[0]} -select-type INDEL -O '$name'.GenotypeGVCFs.raw.INDELs.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantFiltration -R '$genoma' -V '$name'.GenotypeGVCFs.raw.INDELs.vcf.gz -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' -O '$name'.GenotypeGVCFs.hard_filtered.INDELs.vcf.gz
    java -jar /media/storage2/software/Picard/picard.jar MergeVcfs I= '$name'.GenotypeGVCFs.hard_filtered.INDELs.vcf.gz I= '$name'.GenotypeGVCFs.hard_filtered.SNPs.vcf.gz O= '$name'.GenotypeGVCFs.hard_filtered.all.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R '$genoma' -V '$name'.GenotypeGVCFs.hard_filtered.all.vcf.gz --exclude-filtered -O '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz
    vcf-query -l '$name'.GenotypeGVCFs.hard_filtered.fil_rm.all.vcf.gz > GVCF_sample_list       
    """
}
process VCF_to_fasta {
    maxForks 10    
    publishDir "${params.outdir}/snpEff", mode:'copy', pattern: '*.counts.txt'     
    // actualiza una base de datos preexistente
        
    input:    
    each path(GVCF)    
    val(sample_name)
    each path(genome) 
    each path(index_files)
    each path(dict)   

    output:
    path "*.fa", emit: consensus_fasta    
    path "*.counts.txt"
    path "${sample_name}", emit: multiqc_input
        
    script:
    """
    vcftools --gzvcf ${GVCF[0]} --indv ${sample_name} --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > ${sample_name}.vcf 
    cat ${genome} | vcf-consensus ${GVCF[0]} -s ${sample_name}|sed 's/${params.ref_seq_header}/${sample_name}/g' | sed 's/*//g' > ${sample_name}.fa 
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantsToTable -V ${sample_name}.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -F NSAMPLES  -F NCALLED -O ${sample_name}.vcf.sumary_table.tsv 
    grep 'SNP' ${sample_name}.vcf.sumary_table.tsv | wc -l > ${sample_name}.snp.counts.txt
    grep 'INDEL' ${sample_name}.vcf.sumary_table.tsv | wc -l > ${sample_name}.indel.counts.txt       
    sed -i 's/NC_010168.1/Chromosome/g' ${sample_name}.vcf
    java -jar /media/storage2/software/snpEff/snpEff.jar -q -csvStats ${sample_name} -noLog -nodownload Renibacterium_salmoninarum_atcc_33209 ${sample_name}.vcf | vcf-sort > ${sample_name}_snpEff_filtered_snps.vcf    
    """
}
process VCF_to_fasta_GATK_only {
    maxForks 10    
    publishDir "${params.outdir}/snpEff", mode:'copy', pattern: '*.counts.txt'     
    // actualiza una base de datos preexistente
        
    input:    
    each path(GVCF)    
    val(sample_name)
    each path(genome) 
    each path(index_files)
    each path(dict)   

    output:
    path "*.fa", emit: consensus_fasta    
    path "*.counts.txt"
    path "${sample_name}", emit: multiqc_input
        
    script:
    """
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SelectVariants -R ${genome} -V ${GVCF[0]} --sample-name ${sample_name} --exclude-non-variants -O  ${sample_name}.vcf.gz
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar FastaAlternateReferenceMaker -R ${genome} -O ${sample_name}.fa -V ${sample_name}.vcf.gz    
    sed -i 's/^>.*\$/>${sample_name}/g' ${sample_name}.fa
    java -jar /media/storage2/software/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar VariantsToTable -V ${sample_name}.vcf.gz -F CHROM -F POS -F REF -F ALT -F TYPE -F NSAMPLES  -F NCALLED -O ${sample_name}.vcf.sumary_table.tsv 
    grep 'SNP' ${sample_name}.vcf.sumary_table.tsv | wc -l > ${sample_name}.snp.counts.txt
    grep 'INDEL' ${sample_name}.vcf.sumary_table.tsv | wc -l > ${sample_name}.indel.counts.txt
    gunzip -k ${sample_name}.vcf.gz
    sed -i 's/NC_010168.1/Chromosome/g' ${sample_name}.vcf
    java -jar /media/storage2/software/snpEff/snpEff.jar -q -csvStats ${sample_name} -noLog -nodownload Renibacterium_salmoninarum_atcc_33209 ${sample_name}.vcf | vcf-sort > ${sample_name}_snpEff_filtered_snps.vcf   
    """
}
process PhaMe {    
    cpus 30    
    conda '/media/storage2/software/anaconda3/envs/phame_env'
    publishDir "${params.outdir}/PhaMe", mode:'copy'    
    // actualiza una base de datos preexistente
        
    input:    
    path(fasta_files)   
    path(genoma)    
    path(conf_file)    

    output:
    path "workdir/results/trees/", emit: consensus_fasta    
        
    script:
    """  
    mkdir refdir
    cat ${genoma} | sed 's/${params.ref_seq_header}/${params.ref_seq_new_header}/' > ${params.ref_seq_new_header}.fa 
    mv ${params.ref_seq_new_header}.fa refdir/
    mv $fasta_files refdir/    
    mkdir workdir
    /media/storage2/software/PhaME/src/phame $conf_file
    """
}
process Kraken {    
    cpus 30    
    conda '/media/storage2/software/anaconda3/envs/Kraken'  
    maxForks 1  
    tag "Kraken on $sample_id"
    // actualiza una base de datos preexistente
        
    input:    
    tuple val(sample_id), path(reads)  
    
    output:
    path "${sample_id}.report", emit: Kraken_reports    
        
    script:
    """  
    kraken2 --db /media/storage/Adolfo/otros/Trabajos_extras/beagle_UFRJ/metagenomica/output_data/k2_standart_db/ --threads $task.cpus --paired $reads --output ${sample_id}.kraken --report ${sample_id}.report 
    rm ${sample_id}.kraken
    """
}
process Kraken_Biom {    
    cpus 30    
    conda '/media/storage2/software/anaconda3/envs/Kraken'
    // actualiza una base de datos preexistente
        
    input:    
    path(Kraken_reports)  
    
    output:
    path "analisis.biom", emit: Kraken_biom   
        
    script:
    """  
    kraken-biom $Kraken_reports --fmt json -o analisis.biom
    """
}
process Phyloseq {    
    cpus 30    
    conda '/media/storage2/software/anaconda3/envs/phyloseq'
    publishDir "${params.outdir}/Phyloseq", mode:'move'
    // actualiza una base de datos preexistente
        
    input:    
    path(Kraken_reports)  
    
    output:
    path "Analisis_rel_plot.tiff"  
        
    script:
    """
    Rscript /media/storage/Adolfo/otros/Trabajos_extras/renibacteriumSalmoninarum/Phyloseq.R    
    """
}
workflow QC1{
    FASTQC(read_pairs_ch) 
    MULTIQC(FASTQC.out.multiqc_input.collect(),"QC1")
}
workflow Trimm_QC2{
    FASTP(read_pairs_ch)    
    FASTQC(FASTP.out.trimmed_reads)    
    MULTIQC(FASTQC.out.multiqc_input.collect(),"QC2")    
}
workflow index_genome {   
    REF_INDEX(ref_genome_ch)
}
workflow BWA_PhaMe {  
    Trimmed_ch = channel.fromFilePairs(params.Trimmed, checkIfExists: true)
    REF_INDEX(ref_genome_ch)
    BWA_MEM(Trimmed_ch,REF_INDEX.out.indexs,ref_genome_ch,REF_INDEX.out.dict)    
    if ( params.update_joint_vcf == true){
        previus_vcf_ch = channel.fromPath(params.previus_vcf)
        DBImport_update_prev_vcf(BWA_MEM.out.raw_g_vcf.collect(),REF_INDEX.out.indexs,ref_genome_ch,REF_INDEX.out.dict,previus_vcf_ch.collect())        
        GenotypeGVCFs_raw_ch = DBImport_update_prev_vcf.out.GenotypeGVCFs_raw
    } 
    else {
        DBImport(BWA_MEM.out.raw_g_vcf.collect(),REF_INDEX.out.indexs,ref_genome_ch,REF_INDEX.out.dict)
        GenotypeGVCFs_raw_ch = DBImport.out.GenotypeGVCFs_raw
    }
    VCF_filters_GATK_only(GenotypeGVCFs_raw_ch,ref_genome_ch,REF_INDEX.out.indexs,REF_INDEX.out.dict) 
    GVCF_sample_list_ch = VCF_filters_GATK_only.out.sample_list.flatMap{ it.readLines() }
    VCF_to_fasta_GATK_only(VCF_filters_GATK_only.out.vcf,GVCF_sample_list_ch,ref_genome_ch,REF_INDEX.out.indexs,REF_INDEX.out.dict)
    MULTIQC(VCF_to_fasta_GATK_only.out.multiqc_input.collect(),"snpEff")    
    phame_conf_ch = channel.fromPath(params.phame_conf_file)
    PhaMe(VCF_to_fasta_GATK_only.out.consensus_fasta.collect(),ref_genome_ch,phame_conf_ch)
}
workflow Kraken_verification {   
    Trimmed_ch = channel.fromFilePairs(params.Trimmed, checkIfExists: true)
    Kraken(Trimmed_ch)
    Kraken_Biom(Kraken.out.Kraken_reports.collect())
    Phyloseq(Kraken_Biom.out.Kraken_biom)
}
workflow {   
    QC1()
    Trimm_QC2()
    index_genome()
}
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}