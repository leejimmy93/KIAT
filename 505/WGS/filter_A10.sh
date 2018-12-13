#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Julin_share

chrom=`ls *.gz | sed 's/.vcf.gz//g'`
outdir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Julin_share

for i in $chrom
	do 
	echo $i
	vcftools --gzvcf ${i}.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --minQ 40 --minGQ 30 --min-meanDP 1 --max-meanDP 100 --max-missing 0.5 --recode --recode-INFO-all --out ${outdir}/${i}_filtered
	mv ${outdir}/${i}_filtered.recode.vcf ${outdir}/${i}_filtered.vcf 
	gzip ${outdir}/${i}_filtered.vcf 
done
 
