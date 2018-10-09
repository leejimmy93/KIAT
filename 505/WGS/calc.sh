#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_late_silique_no_pop/combined
file_in=$1
prefix=$2

vcftools --gzvcf $file_in \
	 --out $prefix \
         --freq  
echo "done frequency calculation"


vcftools --gzvcf $file_in \
         --out $prefix \
         --depth
echo "done depth calculation"

vcftools --gzvcf $file_in \
         --out $prefix \
         --site-depth
echo "done site depth"

vcftools --gzvcf $file_in \
        --out $prefix \
        --site-mean-depth
echo "done site mean depth"

vcftools --gzvcf $file_in \
         --out $prefix \
         --counts
echo "done counts"

vcftools --gzvcf $file_in \
         --out $prefix \
         --site-pi
echo "done site-pi"

vcftools --gzvcf $file_in \
         --out $prefix \
         --het
echo "done het"

vcftools --gzvcf $file_in \
         --out $prefix \
         --site-quality
echo "done site quality"

vcftools --gzvcf $file_in \
         --out $prefix \
         --missing-indv
echo "done missing indv"

vcftools --gzvcf $file_in \
         --out $prefix \
         --missing-site
echo "done missing site"

vcftools --gzvcf $file_in \
         --out $prefix \
         --hap-r2 \
	 --ld-window-bp 50000
echo "done LD"



