## Adding Samples

There are a few steps necessary to add a sample to the vcf. The script preprocess_bam_to_vcf.sh contains a few commands that might need further testing and modifications but should represent a good starting point.

## Pre-processig

Firstly we need to extract from the big vcf the information about indels within the regions of interest. This can be done using for example bcftools but any other software that deals with this format will do. Here I have used an extended region file (spanning 1kb per region) which is provided in the Input folder.

#### Command
bcftools view -v indels -R chr6_extended_target_list.bed -Ov -o Input/chr6_target_indels.vcf path/to/big/vcf/from/NIH

## Creating Summaries

#### Example Command

## Plotting Results
