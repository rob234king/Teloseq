# TeloseqV2.0: Basic pipeline to analyse telomere enriched data.

HG002 custom reference used and pipeline requires all chr arms to be cut and orientated to 5'->3'.

Note: Future development is for a reference free approach because when you move away from a sample that matches your reference, then you will get mismapping and poor alignments with softclipping due to a repeating unit next to telomere and similarity of some of chr arms. When deviate from HG002 for instance with colo829 then the snps used to separate to chr arm are replaced with others but there is no way of knowing then where these chr arm reads are to be placed by the mapper. So to account for this then a custom reference would be a better option however some arms may not be able to be named using current public reference similarity.

## Prerequisites: 

Nextflow, Java

## Example command on linux system

#no barcodes. Recommended to demultiplex using dorado during basecalling.


nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --input ./test_data/telomerev1.fastq.gz --publishDir test_output -profile conda,mamba

#with barcode demultiplexing if basecalling using bonito custom model. Format is comma separated barcodes, barcode scripts in bin to be expanded as more barcodes added. The number after : is the number of mismatches used with edlib.


nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --input ./test_data/Barcode_test.fastq.gz --publishDir test_output2s --barcode Barcode01:5,Barcode02:5,Barcode05:5 -profile conda,mamba

## Using conda

The pipeline is prepared to manage its software requirements automatically using
conda. If you have it available you can enable this functionality with the
`-profile conda` argument.

### Using mamba

If you also have mamba installed, you can take advantage of its performance
improvements by using `-profile conda,mamba`

## Running in Epi2me labs via Windows on a laptop (next release)

download epi2me labs and go to workflows->import workflows and enter the URL for this github to download the workflow. Docker image to be updated in next release.
https://labs.epi2me.io/downloads/
