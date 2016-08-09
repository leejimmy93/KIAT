#!/bin/bash

directory=`ls -d */ | grep "fastqc"`
echo $directory

for i in $directory
     do

# cat summary file 

cd $i 
cat summary.txt >> /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/fastqc_result/summary_all.txt  

python /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/fastqc_result/extract_inf.py -i fastqc_data.txt -o overrepresented_seq

cat overrepresented_seq >> /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/fastqc_result/overrepresented_seq

cd ..

cat overrepresented_seq | grep -v "#" | awk '{print $1, $4}' | sort | uniq -c | nl | awk '{print ">" $1 "_" $2 "_" $4 "\n" $3}' > overrepresented_seq.fa 

done
