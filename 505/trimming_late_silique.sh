#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/data/late_silique

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i

# trimming 

java -jar /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE Sample_${i}/${i}_1.fq.gz Sample_${i}/${i}_2.fq.gz Sample_${i}/${i}_paired_1.fq.gz Sample_${i}/${i}_unpaired_1.fq.gz Sample_${i}/${i}_paired_2.fq.gz Sample_${i}/${i}_unpaired_2.fq.gz ILLUMINACLIP:Bradseq_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 6

done  

