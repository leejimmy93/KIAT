#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.rapa

mkdir star_genome_V2.1

STAR --runMode genomeGenerate --genomeDir star_genome_V2.1/ --genomeFastaFiles BrapaV2.1PacBio.Chr.fa --sjdbGTFfile BrapaV2.1PacBio.Chr.gene.gff --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS limitGenomeGenerateRAM 61617847680
