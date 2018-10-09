#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw

file=`ls *.gz | grep -v "A01" | sed 's/.gz//g'`

for i in $file	
	do
	echo $i
	vcf-shuffle-cols -t chrA01_random.vcf.gz $i.gz > reordered/$i
	bgzip reordered/$i
	tabix -p vcf reordered/$i.gz
done

cp chrA01_random.vcf.gz reordered
