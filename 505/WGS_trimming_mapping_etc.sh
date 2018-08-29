#!/bin/bash 

# trimming: K314 SJ_142  
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/raw

sample=`echo K314_S441_L003 SJ_142_S447_L003`
echo $sample
trimmed_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/trimmed

for i in $sample
	do
	echo $i
	java -jar /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz $trimmed_dir/${i}/${i}_paired_1.fq.gz $trimmed_dir/${i}/${i}_unpaired_1.fq.gz $trimmed_dir/${i}/${i}_paired_2.fq.gz $trimmed_dir/${i}/${i}_unpaired_2.fq.gz ILLUMINACLIP::2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:76

done

# mapping and convert to bam: K314, SJ_142, SJ_198 
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/trimmed

sample=`echo K314_S441_L003 SJ_142_S447_L003 SJ_198_S480_L003`
mapping_dir=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/mapping

for i in $sample
        do
        echo $i
	bwa mem -t 4 /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa ${i}/${i}_paired_1.fq.gz  ${i}/${i}_paired_2.fq.gz > $mapping_dir/${i}_tmp1.sam
	samtools view $mapping_dir/${i}_tmp1.sam -bT /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus_v4.1.chromosomes.fa -@ 4 > $mapping_dir/${i}_tmp1.bam

done 

# sort: K314, SJ_142, SJ_198, SJ_197, SJ_205, SJ_207 
sample=`echo K314_S441_L003 SJ_142_S447_L003 SJ_198_S480_L003 SJ_197_S481_L003 SJ_205_S489_L003 SJ_207_S487_L003`

for i in $sample
        do
        echo $i
	samtools sort -@ 4 $mapping_dir/${i}_tmp1.bam > $mapping_dir/${i}.bam	       
	rm $mapping_dir/${i}_tmp1.sam
	rm $mapping_dir/${i}_tmp1.bam
done

