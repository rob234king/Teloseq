TeloseqV1: a Basic pipeline to analyse telomere enriched data.

The reference supplied in test_data was created by aligned HG002 reference to chm13 for naming. The original naming of chr arm subtelomeres conflicted with alignments to chm13 similarity so was corrected based upon this. Subsequently is may be best to use your own chr arm reference if you check and disagree with using chm13 namings.

Note: Future development is for a reference free approach because when you move away from a sample that matches your reference, then you will get mismapping and poor alignments with softclipping due to a repeating unit next to telomere and similarity of some of chr arms. When deviate from HG002 for instance with colo829 then the snps used to separate to chr arm are replaced with others but there is no way of knowing then where these chr arm reads are to be placed by the mapper. So to account for this then a custom reference would be a better option however some arms may not be able to be named using current public reference similarity.

#Prerequisites: 

Nextflow, Java

#Example command on linux system

nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta --input ./test_data/test.fastq --publishDir test_output
