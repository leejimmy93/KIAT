#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check/reference_genome_RNAseq

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	vcftools --gzvcf SNP_result/${i}.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minQ 40 --minGQ 30 --min-meanDP 20 --max-meanDP 500 --recode --recode-INFO-all --out SNP_result/${i}_filtered 

done
 
