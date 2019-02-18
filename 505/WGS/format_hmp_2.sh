#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_late_silique_131_sample/combined

file=`echo 505_RNA_impute.hmp.txt`
for i in $file
	do 	
	sed 's/#//g' $i > hmp/${i}
done

	
