#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ol 

# merge bam files from biological replicates
samtools merge Ol_flowering.bam 2/Aligned.sortedByCoord.out.bam All1_Gae_2/Aligned.sortedByCoord.out.bam All1_Gae_3/Aligned.sortedByCoord.out.bam  

echo "done merge"

# extract uniquely mapped reads  
samtools view Ol_flowering.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa - > Ol_flowering_unique.bam

echo "done unique extraction"
# sort bam file 
samtools sort Ol_flowering_unique.bam -o Ol_flowering_unique_sorted.bam

echo "done sort"  
# remove PCR duplicate
samtools rmdup -s Ol_flowering_unique_sorted.bam Ol_flowering_unique_sorted_rmdup.bam 

echo "done remove PCR duplicate" &&
# index bam file 
samtools index Ol_flowering_unique_sorted_rmdup.bam  

echo "done index"
# add read group & call SNPs & include Da-Ae & Da-Ol
bamaddrg -b Ol_flowering_unique_sorted_rmdup.bam -s Ol_flowering_combined -r Ol_flowering_combined | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels > Ol_flowering_combined.vcf

echo "done SNP calling"

 
