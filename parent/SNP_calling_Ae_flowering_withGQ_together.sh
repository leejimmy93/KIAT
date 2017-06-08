#!/bin/bash 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ae
# add read group & call SNPs & include Da-Ae & Da-Ol
bamaddrg -b Ae_flowering_unique_sorted_rmdup.bam -s combined -r combined -b 6_paired.star.trim.dir_unique_sorted_rmdup.bam -s rep1 -r rep1 -b Ae_Gae_2_paired.star.trim.dir_unique_sorted_rmdup.bam -s rep2 -r rep2 -b Ae_Gae_3_paired.star.trim.dir_unique_sorted_rmdup.bam -s rep3 -r rep3 | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa --stdin --dont-left-align-indels --genotype-qualities > Ae_flowering_withGQ_together.vcf


 
