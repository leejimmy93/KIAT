#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/mapping

sample=`echo K314_S441_L003 SJ_142_S447_L003 SJ_198_S480_L003`

for i in $sample
	do
	mkdir $i
	samtools view -F 4 ${i}.bam | grep -v "XA:Z" | grep -v "SA:Z" > ${i}/Unique.sorted.sam
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa ${i}/Unique.sorted.sam > ${i}/Unique.sorted.bam
	samtools rmdup -s ${i}/Unique.sorted.bam ${i}/rmdup.bam
	samtools index ${i}/rmdup.bam
	bamaddrg -b ${i}/rmdup.bam -s ${i} -r ${i} > ${i}/final.bam
	rm ${i}/Unique.sorted.sam
	rm ${i}/Unique.sorted.bam
	rm ${i}/rmdup.bam	
done
