#!/bin/bash

echo Collecing names of all samples

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/B_subgenome_check/output/mapping_result/ABC/bam
files=`ls *.stats`
#echo $files
echo `wc -w <<< $files` Directories found

awk 'BEGIN { print "Sample\tNumber_Unique_Mapped" }' > Stats.tab

for i in $files
    do
    echo Extracting mapping stats from stats file in ${i}
    if [ -f ${i} ]
        then
            echo ${i} > BWA.temp.0
            cat ${i} | head -1 | awk '{print $1}' > BWA.temp.1
            paste BWA.temp.0 BWA.temp.1 > BWA.temp
            cat BWA.temp >> Stats.tab
    else
        echo Error, stats file not found in ${i}
    fi
done

rm BWA.temp*
echo Stats.tab has been created
