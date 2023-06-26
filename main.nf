#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the refence genome
 * and read file can be provided as command line options
 * export PATH="location/miniconda3/bin:$PATH"
 * export PATH="location/java11/bin:$PATH"
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
  //check and unzip reference if needed.
  ref = file(params.reference, checkIfExists: true)
  checkReference(ref)
  uncompressedRefFile = checkReference.out.ref1

  //get stats on telomere reads and filter
  telomere_reads(file_ch)
  //map telomere reads to genome
  mapping(telomere_reads.out.fastq, uncompressedRefFile)
  //filter bam
  filtering(mapping.output)
  //get final telomere arms stats
  results(filtering.out.finalbam, uncompressedRefFile)
  //getVersionInfo()
}


  process checkReference {
    input:
    path ref

    output:
    path "reference.fasta", emit: ref1
    
    script:

    """
    if [[ $ref == *.gz ]]; then
      zcat $ref > reference.fasta
    else
      ln -s $ref reference.fasta
    fi
    """
  }

process telomere_reads {

    conda 'seqtk seqkit'

    input:
    path "input.fastq.gz"

    output:
    path "telomere_filtered.fastq.gz", emit: fastq
    path "raw_data_stats.txt", emit: stats1
    path "Reads_telomere_stats.txt", emit: stats2
    path "Reads_telomere_nonTelomere_stats.txt", emit: stats3
    path "reads_with_cutsites.stats.txt", emit: stats4
    path "telomere_filtered_stats.txt", emit: stats5
    
    publishDir "${params.publishDir}", mode: 'copy', overwrite: false

    script:
    """
    seqkit stats -a input.fastq.gz > raw_data_stats.txt
    seqkit grep -s -R 1:1000 -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" input.fastq.gz > telomere.fastq
    seqkit stats -a telomere.fastq > Reads_telomere_stats.txt
    seqkit grep -v -s -R -1500:-1 -p "TAACCCTAACCCTAACCCTAACCCTAACCC" telomere.fastq > telomereandnon.fastq
    seqkit stats -a telomereandnon.fastq > Reads_telomere_nonTelomere_stats.txt
    seqkit grep -v -s -p "CATCCTCCATCCTC","TACCTCATACCT","CTGAGTCCTGAGTC" telomereandnon.fastq > telomere_filtered.fastq
    seqkit stats -a telomere_filtered.fastq > telomere_filtered_stats.txt
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
    path ref2

    output:
    path "Coverage.csv"
    path "Boxplot_of_Telomere_length.pdf"

    publishDir "${params.publishDir}", mode: 'copy', overwrite: false

    script:
    """
    samtools index $input
    samtools fastq $input | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1" "\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt
    samtools faidx $ref2
    seqkit bam $input 2>filtered.bam.stats
    python ${baseDir}/bin/telomerelengthcov.py filtered.bam.stats locationstelomerelast2.txt ${ref2}.fai

    """
}

/* process getVersionInfo {

    output:
    path 'version_info.txt'
   
    script:
    """
    echo 'Version: ' + getVersion() 
    echo 'Platform: ' + getParams() 
    """
} */
