#!/bin/bash
# Example:  bash extract_inf.sh input_file out_file
file_in=$1
file_out=$2
Label=">>Ov"
while read line
do
	if [ "${line:0:4}" = "$Label" ];then
		while read line
		do
			if [ "${line:0:1}" != ">" ];then
				echo $line >> $file_out
			else
				break
			fi
		done
	fi
done < $file_in
