#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_late_silique_131_sample/combined

# WGS 
# cd WGS_filtered_2

file=`echo 505_RNA_impute.vcf.gz | sed 's/.vcf.gz//g'`
for i in $file	
	do
	echo $i
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -SortGenotypeFilePlugin -inputFile ${i}.vcf.gz -outputFile ${i} -fileType VCF
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -fork1 -vcf ${i}.vcf -export -exportType Hapmap -runfork1
done

