#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/leaf_VS_late_silique

# extract common lines between leaf & late silique 
vcftools --keep leaf_common_ID_92 --gzvcf ../vcf_leaf_no_pop/combined/505.vcf.gz --recode --recode-INFO-all --out 505_common_leaf_92

vcftools --keep late_silique_common_ID_92 --gzvcf ../vcf_late_silique_131_sample/combined/505.vcf.gz --recode --recode-INFO-all --out 505_common_late_silique_92

echo "done extract common individuals"

# filter based on other criteria 
vcftools --gzvcf 505_common_leaf_92.recode.vcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.05 --minQ 40 --min-meanDP 5 --max-meanDP 500 --max-missing 0.5 --minGQ 30 --recode --recode-INFO-all --out 505_filtered_leaf_92

gzip 505_filtered_leaf_92.recode.vcf 

vcftools --gzvcf 505_common_late_silique_92.recode.vcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.05 --minQ 40 --min-meanDP 5 --max-meanDP 500 --max-missing 0.5 --minGQ 30 --recode --recode-INFO-all --out 505_filtered_late_silique_92

gzip 505_filtered_late_silique_92.recode.vcf

echo "done filter SNPs"

# convert to hmp format 
run_pipeline.pl -SortGenotypeFilePlugin -inputFile 505_filtered_leaf_92.recode.vcf.gz -outputFile 505_filtered_leaf_92_sorted -fileType VCF

run_pipeline.pl -Xmx5g -fork1 -vcf 505_filtered_leaf_92_sorted.vcf -export -exportType Hapmap -runfork1

run_pipeline.pl -SortGenotypeFilePlugin -inputFile 505_filtered_late_silique_92.recode.vcf.gz -outputFile 505_filtered_late_silique_92_sorted -fileType VCF

run_pipeline.pl -Xmx5g -fork1 -vcf 505_filtered_late_silique_92_sorted.vcf -export -exportType Hapmap -runfork1

echo "all finished"
