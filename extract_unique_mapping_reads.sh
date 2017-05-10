#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ol_with_ABC

sample=`cat sample_list`
echo $sample

for i in $sample 
	do
	echo ${i} 
	samtools view ${i}/mismatch0/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > ${i}/mismatch0/${i}_Unique.sorted.sam
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa ${i}/mismatch0/${i}_Unique.sorted.sam > ${i}/mismatch0/${i}_Unique.sorted.bam

	rm *.sam 
done 

