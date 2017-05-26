#!/bin/bash 

# combine fastq 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data

gunzip 6_2.fq.gz
gunzip Ae_Gae_2_2.fq.gz
gunzip Ae_Gae_3_2.fq.gz

cat 6_2.fq Ae_Gae_2_2.fq Ae_Gae_3_2.fq > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae/single_end_mapping/merge_fastq_then_trim_mapping_analyze/Ae_flowering.fq

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae/single_end_mapping/merge_fastq_then_trim_mapping_analyze

# zip
gzip Ae_flowering.fq 

# trimming 
java -jar /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE Ae_flowering.fq.gz Ae_flowering_trimmed.fq.gz ILLUMINACLIP:Bradseq_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 6

# mapping 
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome \
        --readFilesIn Ae_flowering_trimmed.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
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
