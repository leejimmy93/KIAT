#!/bin/bash

SRR=`cat SRR_ID`

for i in $SRR
	do
	echo $i
	fastq-dump $i
done
