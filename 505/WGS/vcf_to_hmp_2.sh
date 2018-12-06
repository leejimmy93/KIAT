#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw/filtered/hmp

# WGS 
# cd WGS_filtered_2

file=`ls *.gz | sed 's/_filtered_het_0.2.recode.vcf.gz//g'`
for i in $file	
	do
	echo $i
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -SortGenotypeFilePlugin -inputFile ${i}_filtered_het_0.2.recode.vcf.gz -outputFile ${i} -fileType VCF
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -fork1 -vcf ${i}.vcf -export -exportType Hapmap -runfork1
done

