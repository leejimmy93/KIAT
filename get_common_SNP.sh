#!/bin/bash 

#sample=`ls | grep "recode"`
#for i in $sample
#	do
#	cat $i | grep -v "#" | awk '{print $1 "\t" $2}' > ${i}_pos
#done

# get common position 
sort 6_filtered.recode.vcf_pos Ae_Gae_2_filtered.recode.vcf_pos | uniq -d > 6_Ae_Gae_2
sort Ae_combined_filtered.recode.vcf_pos Ae_Gae_3_filtered.recode.vcf_pos | uniq -d > combined_Ae_Gae_3

sort 6_Ae_Gae_2 combined_Ae_Gae_3 | uniq -d > common_POS

# filter based on position 
sample=`cat sample_list`
echp $sample
for i in $sample
	do
	vcftools --positions common_POS --vcf ${i}_filtered.recode.vcf --recode --recode-INFO-all --out ${i}_common.vcf
done



