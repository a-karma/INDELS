## Adding Samples

There are a few steps necessary to add a sample to the vcf. The script preprocess_bam_to_vcf.sh contains a few commands that might need further testing and modifications but should represent a good starting point.

## Pre-processig

Firstly we need to extract from the big vcf the information about indels within the regions of interest. This can be done using for example bcftools but any other software that deals with this format will do. Here I have used an extended region file (spanning 1kb per region) which is provided in the Input folder.

#### Command
> bcftools view -v indels -R chr6_extended_target_list.bed -Ov -o Input/chr6_target_indels.vcf path/to/big/vcf/from/NIH

## Creating Summaries
Summary statistics can be obtained employing the custom python script (ins_del_summary.py).

This script requires three input arguments: the subsetted vcf file (chr6_target_indels.vcf), a file with individual names one per line in the same order as they appear in the vcf, an integer representing the minimum size of indels to be included in the analysis. Arguments must be provided at the call in this order. 

The output of the script is a headerless tab separated file. Each row represents a single individual and consists of the following fields: 
- Sample ID
- Number of Missing data
- Number of times the REF allele was observed
- Number of Discarded indels due to size filtering
- Number of Insertion events
- Number of Deletion events
- mean length of Insertion in bp
- Standard Deviation of Insertion size
- mean length of Deletion in bp
- Standard Deviation of Deletion size
- Total net length in bp (Deletions from REF are counted as negative, Insertions are considered positive)

Each .tsv output file will be termed Indels_summary_min followed by the value provided as the third argument. These files need to be manually curated to add a stutus column (wild vs domestics) or a functional category for the breed. Examples of these curetad versions are provided in the Output folder (see for example final_summary_min5.tsv)

#### Example Command

> python3 Scripts/ins_del_summary.py Input/chr6_target_indels.vcf Input/ind_names.txt 5

## Plotting Results
