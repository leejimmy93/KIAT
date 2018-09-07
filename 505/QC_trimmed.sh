#!/bin/bash

start=`date +%s`

## Identify each array run
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/trimmed
sample=`echo K314_S441_L003 SJ_142_S447_L003`

for i in $sample
	do echo $i
	fastqc ${i}/${i}_paired_1.fq.gz -o /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/trimmed
	fastqc ${i}/${i}_paired_2.fq.gz -o /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/trimmed
done

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
