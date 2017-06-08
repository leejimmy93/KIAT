#!/bin/bash 

echo "merge bam file"
samtools merge Ae_flowering_single.bam 6/Aligned.toTranscriptome.out.bam Ae_Gae_2/Aligned.toTranscriptome.out.bam Ae_Gae_3/Aligned.toTranscriptome.out.bam 

echo "sort bam file"
samtools sort Ae_flowering_single.bam -o Ae_flowering_single_sorted.bam

echo "extract unique mapping"
samtools view Ae_flowering_single_sorted.bam | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa > Ae_flowering_single_uniq.bam

echo "remove PCR duplicate"
samtools rmdup -s Ae_flowering_single_uniq.bam Ae_flowering_single_rmdup.bam

echo "index bam file"
samtools index Ae_flowering_single_rmdup.bam 

echo "SNP calling"
bamaddrg -b Ae_flowering_single_rmdup.bam -s Ae_floweirng_single -r Ae_flowering_single | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > Ae_flowering_single.vcf

 

