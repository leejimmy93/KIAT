#!/bin/bash 
sample=`ls *.sam | sed 's/_Unique.sorted.sam//g'`

for i in $sample
	do
	echo $i
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus_plus_Bsubgenome/B_napus_plus_Bsub.fa ${i}_Unique.sorted.sam > ${i}_Unique.sorted.bam	
done
