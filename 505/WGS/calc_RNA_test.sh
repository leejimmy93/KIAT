#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_combine_WGS_RNAseq/RNA_filtered

file_in=$1
prefix=$2

vcftools --gzvcf $file_in \
	 --out $prefix \
         --freq  
echo "done frequency calculation"

