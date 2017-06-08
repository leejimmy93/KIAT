#!/bin/bash 

reference=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa

java -cp /usr/local/stow/GenomeAnalysisTK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
    -R $reference \
    -V Ol_realignedBam.bam.xaa.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xab.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xac.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xad.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xae.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xaf.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xag.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xah.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xai.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xaj.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xak.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xal.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xam.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xan.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xao.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xap.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xaq.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xar.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xas.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xat.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xau.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xav.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xaw.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xax.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xay.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xaz.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xba.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbb.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbc.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbd.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbe.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbf.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbg.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbh.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbi.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbj.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbk.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbl.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbm.intervals.vcf.filtered.vcf \
    -V Ol_realignedBam.bam.xbn.intervals.vcf.filtered.vcf \
    -out Ol.vcf \
    -assumeSorted
