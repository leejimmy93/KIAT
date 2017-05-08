#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae 

# merge bam files from biological replicates
sample=`ls | sed 's/\///g'` 

for i in $sample
	do
	echo $i

# extract uniquely mapped reads  
samtools view ${i}/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa - > ${i}_unique.bam

echo "done unique extraction"
# sort bam file 
samtools sort ${i}_unique.bam -o ${i}_unique_sorted.bam

echo "done sort"  
# remove PCR duplicate
samtools rmdup -s ${i}_unique_sorted.bam ${i}_unique_sorted_rmdup.bam 

echo "done remove PCR duplicate" &&
# index bam file 
samtools index ${i}_unique_sorted_rmdup.bam  

echo "done index"
# add read group & call SNPs & include Da-Ae & Da-Ol
bamaddrg -b ${i}_unique_sorted_rmdup.bam -s ${i} -r ${i} | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels > ${i}.vcf

echo "done SNP calling"

done
