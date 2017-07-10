#!/bin/bash 

data_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/data/bioftp.org/TBD170026_20170202
result_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check/use_505_data_check

cd $data_dir

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	echo "extract unique mapping"
	samtools view Sample_${i}/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa > Sample_${i}/unique.bam 
	echo "sort bam file"
	samtools sort Sample_${i}/unique.bam -o Sample_${i}/unique_sorted.bam
	echo "remove PCR duplicate"
	samtools rmdup -s Sample_${i}/unique_sorted.bam Sample_${i}/unique_sorted_rmdup.bam
	echo "index bam file"
	samtools index Sample_${i}/unique_sorted_rmdup.bam
	echo "SNP calling"
	bamaddrg -b Sample_${i}/unique_sorted_rmdup.bam -s ${i} -r ${i} | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > $result_dir/${i}.vcf

done
