#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/check_low_coverage/bioinformatics/BWA/MAF_0.05

# WGS 
# cd WGS_filtered_2

file=`echo chrC03`
for i in $file	
	do
	echo $i
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -SortGenotypeFilePlugin -inputFile ${i}_filtered_het_0.2.vcf.gz -outputFile ${i}_het_0.2 -fileType VCF
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx50g -fork1 -vcf ${i}_het_0.2.vcf -export -exportType Hapmap -runfork1
done

