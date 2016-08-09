#!/bin/bash

file=`ls *.fq`
echo $file

for i in $file
      do

# trimming 
java -jar /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE $i $i.trimmed.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 	if [[ -d $i.dir ]]
#           then
#              echo OK
#        else 
#              mkdir $i.dir
#        fi

#        cd $i.dir/
       ### copy overrepresented seqs to each directory

#        home=`pwd`
#        echo "#!/bin/bash" > $i.sh
#        echo cd $home     >> $i.sh
### trimming  

### CDS mapping
#      echo kallisto quant --single --plaintext -l 250 -s 50 -t 12 -i /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5.cds.19.kai -o /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/$i.dir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/$i >> $i.sh  

### remove 1st row and replace it with "ID & lib ID"
#        (awk '{print $1 "\t" $4}' abundance.tsv | sed 1d | sed '1 i\ID	'$i'' > $i.abundance.tsv) >> $i.sh
       
#       chmod 755 $i.sh
#       ./$i.sh
#        mv $i.abundance.tsv /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/read_count
#       cd ../ 
      
done
 
#        cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/read_count

### combine all abundance.tsv file into one master file  
#       paste *abundance.tsv | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > reads       
#       paste *abundance.tsv | awk '{print $1}' > ID 
#       paste ID reads > read.count.tsv 

