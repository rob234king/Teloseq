# TeloseqV1: Basic pipeline to analyse telomere enriched data.

The reference supplied in test_data was created by aligned HG002 reference to chm13 for naming. The original naming of chr arm subtelomeres conflicted with alignments to chm13 similarity so was corrected based upon this. Subsequently is may be best to use your own chr arm reference if you check and disagree with using chm13 namings.

Note: Future development is for a reference free approach because when you move away from a sample that matches your reference, then you will get mismapping and poor alignments with softclipping due to a repeating unit next to telomere and similarity of some of chr arms. When deviate from HG002 for instance with colo829 then the snps used to separate to chr arm are replaced with others but there is no way of knowing then where these chr arm reads are to be placed by the mapper. So to account for this then a custom reference would be a better option however some arms may not be able to be named using current public reference similarity.

## Prerequisites: 

Nextflow, Java

## Example command on linux system

nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --input ./test_data/telomerev1.fastq.gz --publishDir test_output

## Using conda

The pipeline is prepared to manage its software requirements automatically using
conda. If you have it available you can enable this functionality with the
`-profile conda` argument.

### Using mamba

If you also have mamba installed, you can take advantage of its performance
improvements by using `-profile conda,mamba`

## Running in Epi2me labs via Windows on a laptop

## download epi2me labs and go to workflows->import workflows and enter the URL for this github to download the workflow. 

https://labs.epi2me.io/downloads/
