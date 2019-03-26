#!/bin/bash 

n=`echo 1 2 3 4`
cd ~/bin/structure_kernel_src

for i in $n
	do
	echo $i
	./structure -K 5 -L 148069 -N 238 -i ~/505/WGS/STRUCTURE/505_WGS_pop.txt -o ~/505/WGS/STRUCTURE/505_WGS_pop_K5_${i}
done
