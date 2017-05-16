#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check

sample=`cat sample_list`
echo $sample

for i in $sample
	do

vcftools --gzvcf ${i}_paired.star.trim.dir.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minQ 40 --min-meanDP 20 --max-meanDP 500 --minGQ 30 --recode --recode-INFO-all --out ${i}_filtered 


done
 
