#!/bin/bash 

dup_files=`cat duplicated_file`

for i in $dup_files
	do
	echo $i
	rm /Network/Servers/avalanche.plb.ucdavis.edu${i}
done
