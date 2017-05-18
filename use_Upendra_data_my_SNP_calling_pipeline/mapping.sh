#!/bin/bash 

dir_data=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_fall/BIS180L/trimming_mapping/split_fq
dir_analysis=`pwd`

cd $dir_data 
cat IMB211_NDP_1_PETIOLE.fq IMB211_NDP_2_PETIOLE.fq > tmp1
cat tmp1 IMB211_NDP_3_PETIOLE.fq tmp1 > IMB211.fq 

cat R500_NDP_1_PETIOLE.fq R500_NDP_2_PETIOLE.fq > tmp2
cat R500_NDP_3_PETIOLE.fq tmp2 > R500.fq 

cd $dir_analysis
mkdir IMB211
mkdir R500

mv $dir_data/IMB211.fq IMB211
mv $dir_data/R500.fq R500 

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.rapa/STAR_genome_dir/ \
        --readFilesIn ${i}/${i}.fq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${i}/ \
        --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.rapa/Brapa_gene_v1.5.gff \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --alignIntronMax 15000 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --runThreadN 8 \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon CDS \
        --outReadsUnmapped Fastx  
done
