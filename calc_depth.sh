#!/bin/bash 

sample=`ls *.bam | sed 's/_Unique.sorted.bam//g'`
echo $sample

for i in $sample
do	
	echo $i
	samtools depth ${i}_Unique.sorted.bam -a > ${i}_depth 
	done

