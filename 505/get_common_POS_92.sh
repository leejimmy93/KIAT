#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/leaf_VS_late_silique

# gunzip 505_filtered_late_silique_92.recode.vcf.gz
# gunzip 505_filtered_leaf_92.recode.vcf.gz

# cat 505_filtered_late_silique_92.recode.vcf | grep -v "#" | awk '{print $1 "\t" $2}' > POS_late_silique_92
# cat 505_filtered_leaf_92.recode.vcf | grep -v "#" | awk '{print $1 "\t" $2}' > POS_leaf_92

# vcftools --positions POS_late_silique_92 --vcf 505_filtered_leaf_92.recode.vcf --recode --recode-INFO-all --out 505_filtered_common_leaf_92.vcf
# vcftools --positions POS_leaf_92 --vcf 505_filtered_late_silique_92.recode.vcf --recode --recode-INFO-all --out 505_filtered_common_late_silique_92.vcf

# index 
bgzip -c 505_filtered_common_leaf_92.vcf.recode.vcf > 505_filtered_common_leaf_92.vcf.recode.vcf.gz
bgzip -c 505_filtered_common_late_silique_92.vcf.recode.vcf > 505_filtered_common_late_silique_92.vcf.recode.vcf.gz
tabix -p vcf 505_filtered_common_leaf_92.vcf.recode.vcf.gz
tabix -p vcf 505_filtered_common_late_silique_92.vcf.recode.vcf.gz
