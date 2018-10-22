#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS_RNA

file_in=$1
prefix=$2

vcftools --gzvcf $file_in \
	 --out $prefix \
         --freq  
echo "done frequency calculation"

vcftools --gzvcf $file_in \
        --out $prefix \
        --site-mean-depth
echo "done site mean depth"

vcftools --gzvcf $file_in \
         --out $prefix \
         --site-quality
echo "done site quality"

vcftools --gzvcf $file_in \
         --out $prefix \
         --missing-site
echo "done missing site"

