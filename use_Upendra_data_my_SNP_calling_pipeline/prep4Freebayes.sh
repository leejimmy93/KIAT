#!/bin/bash 

sample=`cat sample_list`

for i in $sample
	do
	echo $i

	samtools view ${i}/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > ${i}.sam
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.rapa/Brapa_sequence_v1.5.fa ${i}.sam > ${i}.bam
	samtools rmdup -s ${i}.bam ${i}_rmdup.bam
	samtools index ${i}_rmdup.bam	
	rm ${i}.sam	
done

# add read group & SNP calling
# bamaddrg -b Ae_no_bolting_unique_mapped_rmdup.bam -s Ae -r Ae -b Ol_no_bolting_unique_mapped_rmdup.bam -s Ol -r Ol | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels > Ae_Ol.vcf


