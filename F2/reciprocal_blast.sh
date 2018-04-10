#!/bin/bash 

# goal: get reciprocal blast between two input blast result files A & B 
# 1) define variable 
A=`ls /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/linolenic_acid/A08_test_genes_2.fa`
B=`ls /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/linolenic_acid/C03_test_genes_2.fa`

echo "$A"
echo "$B"

# 2) blastn with only one taget seq  
blastn -db /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/linolenic_acid/C03_test_genes_2 -query $A -out $A.blast.out -evalue 1e-6 -outfmt 6 -max_target_seqs 1

blastn -db /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/linolenic_acid/A08_test_genes_2 -query $B -out $B.blast.out -evalue 1e-6 -outfmt 6 -max_target_seqs 1  

# 3) get only the seq_id & target_id columns for both outputs of A & B
awk '{print $1 "\t" $2}' $A.blast.out > $A.modified
awk '{print $2 "\t" $1}' $B.blast.out > $B.modified

wc -l $A.modified
wc -l $B.modified 

# 4) compare rows in modified_A with B, if the same, output the result
sort $A.modified $B.modified | uniq -d > reciprocal.blast

echo "done"

