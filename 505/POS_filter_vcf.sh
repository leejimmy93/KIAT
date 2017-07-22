#!/bin/bash 

vcf=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_late_silique_131_sample/combined/505_filtered_0.2_missing_sorted.vcf
pos=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/genomic_prediction/output/SNP.ID.het.filterd.0.2

cat $vcf | grep "#" > header
cat $vcf | grep -v "#" > content 


