#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/flower_Ol
dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/mapping_comparison_A_B_C_subgenome

cd 2/
echo "working on sample 2 from A+B+C"
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/2_ABC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa 2_ABC_Unique.sorted.sam > $dir/2_ABC_Unique.sorted.bam

cd ..

echo "working on sample All1_Gae_2 A+B+C"
cd All1_Gae_2/
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/All1_Gae_2_ABC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa All1_Gae_2_ABC_Unique.sorted.sam > $dir/All1_Gae_2_ABC_Unique.sorted.bam

cd ..

echo "working on sample All1_Gae_3 A+B+C"
cd All1_Gae_3/
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/All1_Gae_3_ABC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa All1_Gae_3_ABC_Unique.sorted.sam > $dir/All1_Gae_3_ABC_Unique.sorted.bam

cd ..

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data

cd 2_paired.star.trim.dir/
echo "working on sample 2 from A+C"
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/2_AC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa 2_AC_Unique.sorted.sam > $dir/2_AC_Unique.sorted.bam

cd ..

cd All1_Gae_2_paired.star.trim.dir/
echo "working on sample All1_Gae_2 A+C"
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/All1_Gae_2_AC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa All1_Gae_2_AC_Unique.sorted.sam > $dir/All1_Gae_2_AC_Unique.sorted.bam

cd ..

echo "working on sample All1_Gae_3 A+C"
cd All1_Gae_3_paired.star.trim.dir/
samtools view Aligned.sortedByCoord.out.bam | awk '$5 == "255"' > $dir/All1_Gae_3_AC_Unique.sorted.sam
samtools view -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa All1_Gae_3_AC_Unique.sorted.sam > $dir/All1_Gae_3_AC_Unique.sorted.bam

cd ..

