#!/bin/bash 

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw/reordered

vcf-concat chrA01_random.vcf.gz chrA01.vcf.gz chrA02_random.vcf.gz chrA02.vcf.gz chrA03_random.vcf.gz chrA03.vcf.gz chrA04_random.vcf.gz chrA04.vcf.gz chrA05_random.vcf.gz chrA05.vcf.gz chrA06_random.vcf.gz chrA06.vcf.gz chrA07_random.vcf.gz chrA07.vcf.gz chrA08_random.vcf.gz chrA08.vcf.gz chrA09_random.vcf.gz chrA09.vcf.gz chrA10_random.vcf.gz chrA10.vcf.gz chrAnn_random.vcf.gz chrC01_random.vcf.gz chrC01.vcf.gz chrC02_random.vcf.gz chrC02.vcf.gz chrC03_random.vcf.gz chrC03.vcf.gz chrC04_random.vcf.gz chrC04.vcf.gz chrC05_random.vcf.gz chrC05.vcf.gz chrC06_random.vcf.gz chrC06.vcf.gz chrC07_random.vcf.gz chrC07.vcf.gz chrC08_random.vcf.gz chrC08.vcf.gz chrC09_random.vcf.gz chrC09.vcf.gz chrCnn_random.vcf.gz chrUnn_random.vcf.gz | gzip -c > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw/combined/505.vcf.gz

echo "done" 
