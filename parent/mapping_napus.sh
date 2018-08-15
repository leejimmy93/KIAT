#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/input
samples=`cat sample_list_Ae`
sample_dir="/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/calenbadger/assembly_parent/B.napus/De_novo_Assembly/Trimmed_Fastqs"
output_dir="/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/output/mapping_result/AC"

echo $samples
for i in $samples
	do 
	echo $i 
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome \
        --readFilesIn $sample_dir/${i}_1.paired.fq.gz  $sample_dir/${i}_2.paired.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $output_dir/${i} \
	--twopassMode Basic \
	--runThreadN 6 \
	--outReadsUnmapped Fastx \
	--readFilesCommand zcat \
	--outFilterMismatchNmax 0

done
