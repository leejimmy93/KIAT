#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae 

samtools view Ae_Gae_3_paired.star.trim.dir/Aligned.sortedByCoord.out.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa > Ae_Gae_3_paired.star.trim.dir_unique.bam

samtools sort Ae_Gae_3_paired.star.trim.dir_unique.bam -o Ae_Gae_3_paired.star.trim.dir_unique_sorted.bam

samtools rmdup -s Ae_Gae_3_paired.star.trim.dir_unique_sorted.bam Ae_Gae_3_paired.star.trim.dir_unique_sorted_rmdup.bam

samtools index Ae_Gae_3_paired.star.trim.dir_unique_sorted_rmdup.bam 


