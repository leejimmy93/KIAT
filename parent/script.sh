cat Ae_bwa_mem.sam | awk '$3 == "*"{print $1}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/assembly_parent/De_novo/Bwa_B.napus/Ae_napus_unmapped_ID

cat Ol_bwa_mem.sam | awk '$3 == "*"{print $1}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/assembly_parent/De_novo/Bwa_B.napus/Ol_napus_unmapped_ID

# tag uniquely mapped & multiple mapped genes 
cat Ae_bwa_mem.sam | awk 'NF > 8' | grep -v "XA:Z" | awk '{print $1}' | sort | uniq  | awk '{print $0}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/assembly_parent/De_novo/Bwa_B.napus/Ae_napus_unique_mapped_ID

cat Ol_bwa_mem.sam | awk 'NF > 8' | grep -v "XA:Z" | awk '{print $1}' | sort | uniq  | awk '{print $0}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/assembly_parent/De_novo/Bwa_B.napus/Ol_napus_unique_mapped_ID

# check Uniref hit orgamism 
cat Ae_uniprot | awk 'BEGIN{FS="\t"}{print $7}' | sort | uniq -c
cat Ol_uniprot | awk 'BEGIN{FS="\t"}{print $7}' | sort | uniq -c
