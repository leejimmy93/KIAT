#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/output/mapping_result/ABC/bam

sample=`ls | grep ^6 | grep -v "bai" | tail -49`

for i in $sample
	do
	echo $i
	samtools flagstat ${i} > ${i}.stats
done 
