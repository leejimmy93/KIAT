#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS_RNA

file=`ls *.vcf | grep -v "A10"`

for i in $file
	do
	echo $i
	bash /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/KIAT/505/WGS/calc_WGS_test.sh $i $i 
done
