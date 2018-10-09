#!/bin/bash

# compress and index partial vcf files
echo "compress and index partial vcf files"
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw

sample=`ls *.vcf`
for i in $sample
	do
        echo $i
        bgzip $i
        tabix -p vcf $i.gz
done
