#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the refence genome
 * and read file can be provided as command line options
 * export PATH="/home/OXFORDNANOLABS/rking/miniconda3/bin:$PATH"
 * export PATH="/home/OXFORDNANOLABS/rking/miniconda3/envs/java11/bin:$PATH"
 */

// Declare syntax version
nextflow.enable.dsl=2


params.input = ''
params.reference = ''
params.publishDir = "TeloseqV1"

println """\

         O N T   T E L O S E Q    P I P E L I N E  (P R O T O T Y P E : B A S I C v1)
         =================================================================
         Public reference   : ${params.reference}
         Fastq read files   : ${params.input}
         Publish directory  : ${params.publishDir}
         """
         .stripIndent()

file_ch = Channel.fromPath("${params.input}")

workflow {
  telomere_reads(file_ch)
  ref = file(params.reference, checkIfExists: true)
  mapping(telomere_reads.out.fastq, ref)
  filtering(mapping.output)
  results(filtering.out.finalbam, ref)
}

process telomere_reads {

    conda 'seqtk seqkit'

    input:
    file "input.fastq.gz"

    output:
    path "telomere_filtered.fastq.gz", emit: fastq
    path "raw_data_stats.txt", emit: stats1
    path "Reads_telomere_stats.txt", emit: stats2
    path "Reads_telomere_nonTelomere_stats.txt", emit: stats3
    path "reads_with_cutsites.stats.txt", emit: stats4
    

    publishDir "${params.publishDir}", mode: 'copy', overwrite: false

    script:
    """
    seqkit stats input.fastq.gz > raw_data_stats.txt
    seqkit grep -s -R 1:1000 -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" input.fastq.gz > telomere.fastq
    seqkit stats telomere.fastq > Reads_telomere_stats.txt
    seqkit grep -v -s -R -1500:-1 -p "TAACCCTAACCCTAACCCTAACCCTAACCC" telomere.fastq > telomereandnon.fastq
    seqkit stats telomereandnon.fastq > Reads_telomere_nonTelomere_stats.txt
    seqkit grep -v -s -p "CATCCTCCATCCTC","TACCTCATACCT","CTGAGTCCTGAGTC" telomereandnon.fastq > telomere_filtered.fastq
    gzip telomere_filtered.fastq
    seqkit grep -s -p "GATATC" telomere_filtered.fastq.gz | seqkit stats -a - > reads_with_cutsites.stats.txt
    """
}

process mapping {

    conda 'minimap2 samtools'

    input:
    path input
    path ref

    output:
    path "telomere_filter_q10.bam"

    script:
    """
    minimap2 -ax map-ont -t 20 $ref $input | samtools sort -@2 -o telomere_filter.bam
    samtools index telomere_filter.bam
    samtools view -bq 10 -h telomere_filter.bam > telomere_filter_q10.bam
    samtools index telomere_filter_q10.bam
    """
}

process filtering {

    conda 'samtools pandas numpy seqkit'

    input:
    path inputbam

    output:
    path "filtered.bam", emit: finalbam
    path "filtered.bam.bai"

    publishDir "${params.publishDir}", mode: 'copy', overwrite: false

    script:
    """
    samtools index $inputbam
    seqkit bam $inputbam 2>telomereCustomRef_q10.bam.stats
    python ${baseDir}/bin/Identifyreadstoremove.py telomereCustomRef_q10.bam.stats
    samtools view -N ID.txt -U filtered.bam -o /dev/null $inputbam
    samtools index filtered.bam
    """
}


process results {

    conda 'samtools pandas numpy coreutils matplotlib seqkit'

    input:
    path input
    path "reference.fasta"

    output:
    file "Coverage.csv"
    file "Boxplot_of_Telomere_length.pdf"

    publishDir "${params.publishDir}", mode: 'copy', overwrite: false

    script:
    """
    samtools index $input
    samtools fastq $input | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1" "\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt
    samtools faidx reference.fasta
    seqkit bam $input 2>filtered.bam.stats
    python ${baseDir}/bin/telomerelengthcov.py filtered.bam.stats locationstelomerelast2.txt reference.fasta.fai

    """
}


