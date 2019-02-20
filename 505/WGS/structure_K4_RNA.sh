#!/bin/bash 

n=`echo 1 2 3 4`
cd ~/bin/structure_kernel_src

for i in $n
	do
	echo $i
	./structure -K 4 -L 53895 -N 131 -i ~/505/vcf_late_silique_131_sample/combined/STRUCTURE/505_RNA_pop.txt -o ~/505/vcf_late_silique_131_sample/combined/STRUCTURE/505_RNA_pop_K4_i
done
