#!/bin/bash 

samtools sort R500_rmdup.bam -o R500_rmdup.sorted.bam
samtools sort IMB211_rmdup.bam -o IMB211_rumdup.sorted.bam

samtools index R500_rmdup.sorted.bam
samtools index IMB211_rumdup.sorted.bam 
