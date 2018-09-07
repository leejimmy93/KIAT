#!/bin/bash

# module load samtools/1.5
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/mapping

sample=`echo K314_S441_L003 SJ_142_S447_L003 SJ_198_S480_L003`

for i in $sample
	do
	echo $i
	samtools flagstat ${i}.bam > ${i}/${i}.stats
	samtools flagstat ${i}/final.bam > ${i}/final.stats
done
