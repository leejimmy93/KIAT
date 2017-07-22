#!/bin/bash 
# eg ./script $file_in $file_mid
file_in=$1
file_mid=$2

run_pipeline.pl -SortGenotypeFilePlugin -inputFile $file_in -outputFile $file_mid -fileType VCF

run_pipeline.pl -Xmx5g -fork1 -vcf ${file_mid}.vcf -export -exportType Hapmap -runfork1

rm $file_mid


