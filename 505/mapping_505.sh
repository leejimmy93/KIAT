#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/data/late_silique
samples=`cat sample_list`

echo $samples
for i in $samples
	do 
	echo $i 
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome \
        --readFilesIn Sample_${i}/${i}_paired_1.fq.gz  Sample_${i}/${i}_paired_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix Sample_${i}/ \
        --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --alignIntronMax 15000 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --runThreadN 8 \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon CDS \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat

done
