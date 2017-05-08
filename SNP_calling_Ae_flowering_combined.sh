#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae 

# merge bam files from biological replicates
samtools merge Ae_flowering.bam 6_paired.star.trim.dir/Aligned.sortedByCoord.out.bam Ae_Gae_2_paired.star.trim.dir/Aligned.sortedByCoord.out.bam Ae_Gae_3_paired.star.trim.dir/Aligned.sortedByCoord.out.bam  

echo "done merge"

# extract uniquely mapped reads  
samtools view Ae_flowering.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa - > Ae_flowering_unique.bam

echo "done unique extraction"
# sort bam file 
samtools sort Ae_flowering_unique.bam -o Ae_flowering_unique_sorted.bam

echo "done sort"  
# remove PCR duplicate
samtools rmdup -s Ae_flowering_unique_sorted.bam Ae_flowering_unique_sorted_rmdup.bam 

echo "done remove PCR duplicate" &&
# index bam file 
samtools index Ae_flowering_unique_sorted_rmdup.bam  

echo "done index"
# add read group & call SNPs & include Da-Ae & Da-Ol
bamaddrg -b Ae_flowering_unique_sorted_rmdup.bam -s Ae_flowering_combined -r Ae_flowering_combined | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels > Ae_flowering_combined.vcf

echo "done SNP calling"

 
