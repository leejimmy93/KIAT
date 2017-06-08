#!/bin/bash 

sample=`cat sample_list`

echo $sample

for i in $sample
	do
	echo $i
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome \
        --readFilesIn ../../${i}_paired_1.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${i}/ \
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
