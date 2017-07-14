#!/bin/bash 

data_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check/reference_genome_RNAseq
result_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check/reference_genome_RNAseq/SNP_result

cd $data_dir

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	echo "extract unique mapping"
	samtools view ${i}/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa > ${i}/unique.bam 
	echo "sort bam file"
	samtools sort ${i}/unique.bam -o ${i}/unique_sorted.bam
	echo "remove PCR duplicate"
	samtools rmdup -s ${i}/unique_sorted.bam ${i}/unique_sorted_rmdup.bam
	echo "index bam file"
	samtools index ${i}/unique_sorted_rmdup.bam
	echo "SNP calling"
	bamaddrg -b ${i}/unique_sorted_rmdup.bam -s ${i} -r ${i} | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > $result_dir/${i}.vcf
	rm $data_dir/${i}/unique.bam
	rm $data_dir/${i}/unique_sorted.bam
	rm $data_dir/${i}/unique_sorted_rmdup.bam

done
