#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_filtered/two_options

file=`echo 505_filtered_het_0.hmp.txt`
for i in $file
	do 	
	sed 's/#//g' $i > hmp/${i}
done

	