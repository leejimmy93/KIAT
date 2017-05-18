#!/bin/bash 

vcftools --gzvcf R500_IMB211.vcf --remove-indels --min-alleles 2 --max-alleles 2 --minQ 40 --min-meanDP 2 --max-meanDP 500 --max-missing 0.5 --minGQ 30 --recode --recode-INFO-all --out R500_IMB211 

mv R500_IMB211.recode.vcf 505_filtered.vcf
 
