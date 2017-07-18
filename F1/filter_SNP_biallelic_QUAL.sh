#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F1_SNP/vcf/with_GQ

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	vcftools --gzvcf ${i}.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minQ 40 --recode --recode-INFO-all --out ${i}_filtered 

done
 
