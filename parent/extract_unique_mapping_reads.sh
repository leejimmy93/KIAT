#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/output/mapping_result/AC

sample=`cat sample_list_Ae`
echo $sample

for i in $sample 
	do
	echo ${i} 
	samtools view ${i}Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > ${i}_Unique.sorted.sam
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa ${i}_Unique.sorted.sam > ${i}_Unique.sorted.bam

	rm *.sam 
done 

