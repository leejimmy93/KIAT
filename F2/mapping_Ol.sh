#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ol
samples=`ls | sed 's/\///g'`

echo $samples
for i in $samples
	do 
	echo $i 
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus_plus_Bsubgenome/star_genome \
        --readFilesIn ${i}/${i}_paired_1.fq.gz  ${i}/${i}_paired_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${i}/ \
        --twopassMode Basic \
        --runThreadN 12 \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat

done
