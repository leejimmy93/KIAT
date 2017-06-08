#!/bin/bash

# the purpose is to check the existence of B subgenome in Da-Ol
# although a previous run had 0 reads mapped to B subgenome, running a slitely different pipeline found all 505 have some B subgenome, so I need to try that pipeline on Da-Ol 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/mapping_comparison_A_B_C_subgenome/mismatch_0/unique_mapped/2/new

sample=Aligned.sortedByCoord.out.bam
chrom=`cat chrom_list_napus_plus_Bsub`

# get unique bam file, index unique bam file & subset by chrom id & calculate per base depth 
for i in $sample
	do
	echo $i
	samtools view $i | awk '$5 == "255"' > ${i}_unique
	samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus_plus_Bsubgenome/B_napus_plus_Bsub.fa ${i}_unique > ${i}_unique.bam		
	samtools index ${i}_unique.bam
		for j in $chrom
		do
		echo $j
        	samtools view -@ 2 -b ${i}_unique.bam ${j} > ${j}.bam
        	samtools depth ${j}.bam -a > ${j}.depth  
	done

rm ${i}_unique
rm ${i}_unique.bam 

done


