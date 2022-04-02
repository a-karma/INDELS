# Preprocessing
Firstly, we need to create two input files for bedtools to calculate coverage.
The first one (termed Can_Fam_bedtools_genome_file.txt) contains the chromosome sizes of our reference genome (canFam3.1).
The second one (termed chr6_target.bed) is a tab separated file in bed format (Browser Extensible Data) which contains info about the regions of interest.

In order to create Can_Fam_bedtools_genome_file.txt we can make use of the header of any bam aligned to the ref genome by running:
### Commands
samtools view -H <an_input_file.bam> | grep '@SQ' | awk 'BEGIN {OFS="\t"};{print $2,$3}' | sed 's/SN://' | sed 's/LN://' > Can_Fam_bedtools_genome_file.txt
