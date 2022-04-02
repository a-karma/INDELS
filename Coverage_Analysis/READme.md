# Coverage Tables
Firstly, we need to create two input files for bedtools to calculate coverage.

The first one (termed Can_Fam_bedtools_genome_file.txt) contains the chromosome sizes of our reference genome (canFam3.1). The second one (termed chr6_target.bed) is a tab separated file in bed format (Browser Extensible Data) which contains info about the regions of interest.

In order to create Can_Fam_bedtools_genome_file.txt we can make use of the header of any bam aligned to the ref genome while chr6_target.bed was created manually and it is provided in the input folder. For each genome included in the analysis, raw coverage tables were obtained using bedtools. 
### Commands
samtools view -H <any_input_file.bam> | grep '@SQ' | awk 'BEGIN {OFS="\t"};{print $2,$3}' | sed 's/SN://' | sed 's/LN://' > ./Input/Can_Fam_bedtools_genome_file.txt

bedtools coverage -a ./Input/chr6_target.bed -b input.bam -hist -sorted -g ./Input/Can_Fam_bedtools_genome_file.txt > ./Input/raw/input.cov.tsv

# Aggregating Tables
Once all tables have been produced, it is possible to aggregate them into a single dataframe using the custom python script provided (extract_aggregate.py). The dataframe was then edited manually to add the status column (wild vs domestic). 
### Commands
ls ./Input/raw/*.tsv > ./Input/tables_list.tsv

python3.8 ./Scripts/extract_aggregate.py ./Input/chr6_target.bed ./Input/table_list.txt

# Analysis
This edited dataframe (edited_coverage_database.csv) serves as input for the heatmap which shows a per region summary of the breadth of coverage (BoC) and the coverage depth (DoC) for each sample. 

The dataframe is then split by region and each file is then processed using the R script provided (
