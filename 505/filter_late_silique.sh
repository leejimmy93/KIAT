#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_late_silique_131_sample/combined/

vcftools --gzvcf 505.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.05 --minQ 40 --minGQ 30 --min-meanDP 5 --max-meanDP 500 --max-missing 0.2 --recode --recode-INFO-all --out 505_filtered 

mv 505_filtered.recode.vcf 505_filtered_0.2_missing.vcf
gzip 505_filtered_0.2_missing.vcf
 
