#############
# GATK test #
#############

# Step zero: downloading, cleaning and indexing bams
wget -i https://tigress-web.princeton.edu/~vonholdt/WBS_bam_phased/
rm compressed*
rm *.gif*
rm index.html*
ls *.gz > bams
while read -r line; do gunzip $line; done<bams
rm bams
ls *.bam > bam_list.txt
while read -r line; do samtools index $line; done<bam_list.txt
conda activate bamtools
cat Bridgets_bams/bam_list.txt | xargs -P 5 -L 1 ./scripts/add_RG_tag.sh
conda deactivate
cd Bridgets_bams; ls * RG_tag.bam > new_list.txt 
while read -r line; do samtools index $line; done<new_list.txt

# Step 1: running Haplotype Caller
cd /data/shared/REF/CANIS/
gatk CreateSequenceDictionary -R Canis_familiaris.CanFam3.1.dna.toplevel.fa
cd /data/alberto/Will_Synd/
mkdir gvcf
cat Bridgets_bams/new_list.txt | xargs -P 5 -L1 ./scripts/haplotypecaller.sh

# Step 2: Variant Calling
# useful cmd for inspection 
bcftools view --types indels 2782_RG_tag.bam.g.vcf.gz

# Will have then to merge those individual gvcf with the big vcf (bcftools should suffice)
