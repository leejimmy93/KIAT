#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/vcf_combine_WGS_RNAseq

# WGS 
cd WGS_filtered

file=`ls *.gz | sed 's/_filtered.vcf.gz//g'`
for i in $file	
	do
	echo $i
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile ${i}_filtered.vcf.gz -outputFile ${i}_filtered -fileType VCF
	/usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx5g -fork1 -vcf ${i}_filtered.vcf -export -exportType Hapmap -runfork1
done

cd ../RNA_filtered

file=`ls *.gz | sed 's/_filtered.vcf.gz//g'`
for i in $file
        do
        echo $i
        /usr/local/stow/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile ${i}_filtered.vcf.gz -outputFile ${i}_filtered -fileType VCF
        /usr/local/stow/tassel-5-standalone/run_pipeline.pl -Xmx5g -fork1 -vcf ${i}_filtered.vcf -export -exportType Hapmap -runfork1
done
