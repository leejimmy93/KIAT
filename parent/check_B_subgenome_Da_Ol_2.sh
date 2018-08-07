#!/bin/bash

# the purpose is to check the existence of B subgenome in Da-Ol
# although a previous run had 0 reads mapped to B subgenome, running a slitely different pipeline found all 505 have some B subgenome, so I need to try that pipeline on Da-Ol 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/output/mapping_result/AC

sample=`ls | grep "Unique"`
chrom=`cat chrom_list`

# index unique bam file & subset by chrom id & calculate per base depth 
for i in $sample
	do
	echo $i
	# samtools view $i | awk '$5 == "255"' > ${i}_unique
	# samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus_plus_Bsubgenome/B_napus_plus_Bsub.fa ${i}_unique > ${i}_unique.bam		
	samtools index $i
		for j in $chrom
		do
		echo $j
        	samtools view -@ 2 -b $i ${j} > ${i}_${j}.bam
        	samtools depth ${i}_${j}.bam -a > ${i}_${j}.depth  
	done

rm ${i}_${j}.bam 

done


