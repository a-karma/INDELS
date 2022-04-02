# Coverage Tables
Firstly, we need to create two input files for bedtools to calculate coverage.

The first one (termed Can_Fam_bedtools_genome_file.txt) contains the chromosome sizes of our reference genome (canFam3.1). The second one (termed chr6_target.bed) is a tab separated file in bed format (Browser Extensible Data) which contains info about the regions of interest.

In order to create Can_Fam_bedtools_genome_file.txt we can make use of the header of any bam aligned to the ref genome while chr6_target.bed was created manually and it is provided in the input folder. For each genome included in the analysis, raw coverage tables were obtained using bedtools. 
### Commands
samtools view -H <any_input_file.bam> | grep '@SQ' | awk 'BEGIN {OFS="\t"};{print $2,$3}' | sed 's/SN://' | sed 's/LN://' > ./Input/Can_Fam_bedtools_genome_file.txt

bedtools coverage -a ./Input/chr6_target.bed -b input.bam -hist -sorted -g ./Input/Can_Fam_bedtools_genome_file.txt > ./Input/raw/input.cov.tsv

# Aggregating Tables
