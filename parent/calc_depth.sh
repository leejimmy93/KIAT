#!/bin/bash 

sample=`ls | grep "bam"`
echo $sample

for i in $sample
do	
	echo $i
	samtools depth ${i} -a > ${i}_depth 
	done

