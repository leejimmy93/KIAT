#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae 

# merge bam files from biological replicates
sample=`cat sample_list` 

for i in $sample
	do
	echo $i

# add read group & call SNPs & include Da-Ae & Da-Ol
bamaddrg -b ${i}_unique_sorted_rmdup.bam -s ${i} -r ${i} | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > ${i}.vcf

echo "done SNP calling"

done
