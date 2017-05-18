#!/bin/bash 

# add read group & SNP calling
bamaddrg -b R500_rmdup.sorted.bam -s R500 -r R500 -b IMB211_rumdup.sorted.bam -s IMB211 -r IMB211 | freebayes_v0.9.21-7-g7dd41db --fasta-reference /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.rapa/Brapa_sequence_v1.5.fa --stdin --dont-left-align-indels  --genotype-qualities > R500_IMB211.vcf


