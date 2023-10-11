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
params.barcode = ''

println """\

         O N T   T E L O S E Q    P I P E L I N E  (P R O T O T Y P E : B A S I C v2)
         =================================================================
         Public reference   : ${params.reference}
         Fastq read files   : ${params.input}
         Publish directory  : ${params.publishDir}
         Barcode:mismatches  : ${params.barcode}
         """
         .stripIndent()

file_ch = Channel.fromPath("${params.input}")
if( params.barcode ){
         barcodeset="${params.barcode}"
         //print(barcodeset)
         /////////////////////////////////////////////
         //get CSV parameter for script
         barcodesMap = [:]
         // Split the input barcodes into individual items
         barcodesList = "${params.barcode}".split(',')
         // // Parse the barcode information and create a map
         barcodesList.each { barcodeInfo ->
             def (barcode, count) = barcodeInfo.split(':')
             barcodesMap[barcode] = count.toInteger()
         }
         // Generate the list of barcode names
         barcodeNames = barcodesMap.keySet()
         //print(outputFileNames)
         //print(barcodeNames)
         barcodeNamesString = barcodeNames.collect { "${it}.csv " }.join()
         print (barcodeNamesString)
         /////////////////////////////////////////////////
}

workflow {
  //check and unzip reference if needed.
  ref = file(params.reference, checkIfExists: true)
  checkReference(ref)
  uncompressedRefFile = checkReference.out.ref1

 //get stats on telomere reads and filter
  telomere_reads(file_ch)

  //barcode switch if barcode parameter given and barcode seperate if so.
  if( params.barcode ){
    barcode_separate(telomere_reads.out.fastq2)
    //flatten to create seperate set of items to be run 
    newchannel2=barcode_separate.out.barcodefastq.flatten()
    mappingbam(newchannel2,uncompressedRefFile)
  }
  else{
    //map telomere reads to genome
    mappingbam(file_ch, uncompressedRefFile)
  }

  //filter bam
  filtering(mappingbam.output, uncompressedRefFile, telomere_reads.out.baseerror1)

  //get final telomere arms stats
  results(filtering.out.finalbam, uncompressedRefFile)

  //tool versions
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

    conda 'seqtk seqkit samtools pandas'

    input:
    path "input.fastq.gz"

    output:
    path "telomere_filtered.fastq.gz", emit: fastq
    path "telomere_filtered.fastq", emit: fastq2
    path "raw_data_stats.txt"
    path "Basecalling_error_motif_5_500.txt", emit: baseerror1
    
    publishDir "${params.publishDir}/Telomere/", mode: 'copy', overwrite: false

    script:
    """
    seqkit stats -a input.fastq.gz > raw_data_stats.txt
    seqkit grep -v -s -R -60:-1 -p "TAACCCTAACCCTAACCCTAACCCTAACCC" input.fastq.gz > telomere_filtered.fastq
    seqkit locate --only-positive-strand -m 0 -p TATAGT,TATACA telomere_filtered.fastq > RK_motif.txt
    python ${baseDir}/bin/Filter_motif_5_within500.py > Basecalling_error_motif_5_500.txt
    gzip -c telomere_filtered.fastq > telomere_filtered.fastq.gz
    """
}


process barcode_separate {

    conda 'seqkit pandas numpy pyfastx edlib biopython'

    input:
    path input

    output:
    path "*exclusive.txt.fastq", emit: barcodefastq
    path "*csv", emit: barcodecsv

    publishDir "${params.publishDir}/Barcoded/", mode: 'copy', overwrite: false

    script:
    """
    python ${baseDir}/bin/Barcode_identificationv2.py $input $barcodeset
    python ${baseDir}/bin/barcode_IDV2.py $barcodeNamesString
    for i in unique_names_*_exclusive.txt; do
    seqkit grep -f \$i telomere_filtered.fastq > \$i.fastq;done

    """
}




process mappingbam {

    conda 'minimap2 samtools'

    input:
    path fastqFile
    path ref

    output:
    path "${fastqFile.simpleName}.bam"
    publishDir "${params.publishDir}/Raw_bam/", mode: 'copy', overwrite: false

    script:
    """

    minimap2 -ax map-ont -t 20 $ref $fastqFile | samtools sort -@2 -o ${fastqFile.simpleName}.telomere_filter.bam
    samtools index ${fastqFile.simpleName}.telomere_filter.bam
    samtools view -bq 10 -h ${fastqFile.simpleName}.telomere_filter.bam > ${fastqFile.simpleName}.bam
    samtools index ${fastqFile.simpleName}.bam
    """
}

process filtering {

    conda 'samtools pandas numpy seqkit pysam'

    input:
    path inputbam
    path ref2
    path baseerror

    output:
    path "${inputbam.baseName}.filtered.bam", emit: finalbam
    path "${inputbam.baseName}.filtered.bam.bai"

    publishDir "${params.publishDir}/Mapped_telomere/", mode: 'copy', overwrite: false

    script:
    """
    samtools index $inputbam
    seqkit locate --only-positive-strand -m 0 -p GATATC $ref2 > cutsites2.txt
    awk '!a[\$1]++' cutsites2.txt > cutsites3.txt
    awk -F'\t' '{print \$1"\t"\$5"\t"\$5}' cutsites3.txt > cutfirst5.txt
    grep -v 'seqID' cutfirst5.txt > cutfirst7.bed

    seqkit bam $inputbam 2>${inputbam.baseName}.txt

    seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC $ref2 > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1"\t"\$5"\t"\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt
    grep -v 'seqID' locationstelomerelast2.txt > locationstelomerelast2.bed
    
    python ${baseDir}/bin/Read_IDs_to_keep.py ${inputbam.baseName}.txt $baseerror cutfirst7.bed locationstelomerelast2.bed id.txt

    awk '!seen[\$0]++' id.txt > id2.txt 
    samtools view -N id2.txt -o ${inputbam.baseName}.temp.bam $inputbam 
    samtools index ${inputbam.baseName}.temp.bam
    python ${baseDir}/bin/duplicate_removal.py ${inputbam.baseName}.temp.bam ${inputbam.baseName}.dedup.bam
    samtools index ${inputbam.baseName}.dedup.bam
    samtools view -hq 5 -o ${inputbam.baseName}.filtered.bam ${inputbam.baseName}.dedup.bam
    samtools index ${inputbam.baseName}.filtered.bam

    """
}


process results {

    conda 'pandas numpy coreutils matplotlib samtools seqkit'

    input:
    path input
    path ref2

    output:
    path "${input.baseName}_Coverage.csv"
    path "${input.baseName}_Boxplot_of_Telomere_length.pdf"
    path "${input.baseName}_Per_Read_telomere_length.csv"

    publishDir "${params.publishDir}/Summary_Results", mode: 'copy', overwrite: false

    script:
    """
    samtools index $input
    samtools fastq $input | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1" "\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt
    samtools faidx $ref2
    seqkit bam $input 2>${input.baseName}
    python ${baseDir}/bin/telomerelengthcov.py ${input.baseName} locationstelomerelast2.txt ${ref2}.fai

    """
}

// process getVersionInfo {

//     output:
//     path 'version_info.txt'

//     publishDir "${params.publishDir}/Report", mode: 'copy', overwrite: false
//     script:
//     """
//     echo 'Version: ' + getVersion()  > version_info.txt
//     echo 'Platform: ' + getParams() >> version_info.txt
//     """
// } 
