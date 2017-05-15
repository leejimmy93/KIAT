#!/bin/bash 

sample=`cat sample_list`
echo $sample

for i in $sample
	do
	echo $i
	sed 's/twopts.LOD3_rf0.5.Rdata/twopts.f2.05.12.Rdata/g' $i > ${i}_tmp
	sed 's/twopts.f2.LOD3_rf0.5/twopts.f2.05.12/g' ${i}_tmp > ${i}_final	
done
