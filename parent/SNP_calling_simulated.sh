#!/bin/bash 

input=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/jmaloof/2015/rnaseq_simulations/Bnapus_sim/STAR_Bnapus_sim1Aligned.out.SORT.bam
result_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Da-Ae_Da_Ol-1_heterzygosity_check/Julin_simulation

echo "extract unique mapping"
samtools view ${input} | awk '$5 == "255"' | samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa > $result_dir/unique.bam

echo "sort bam file"
samtools sort $result_dir/unique.bam -o $result_dir/unique_sorted.bam

echo "remove PCR duplicate"
samtools rmdup -s $result_dir/unique_sorted.bam $result_dir/unique_sorted_rmdup.bam

echo "index bam file"
samtools index $result_dir/unique_sorted_rmdup.bam

echo "SNP calling"
bamaddrg -b $result_dir/unique_sorted_rmdup.bam -s simulated -r simulated | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > $result_dir/simulated.vcf
 

 
