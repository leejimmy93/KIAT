#!/bin/bash 

chrom=`cat chrom_list_napus_plus_Bsub`
for i in $chrom 
	do
	echo $i
	grep "$i" 2_ABC_depth > 2_ABC_depth_${i} 
	grep "$i" 2_AC_depth > 2_AC_depth_${i}

done 
