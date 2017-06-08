#!/bin/bash 

file=`ls *.vcf`

for i in $file 
	do 
	GATK -T VariantFiltration -R ~/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa -V $i -window 35 -cluster 3 -o $i.filtered.vcf
done
