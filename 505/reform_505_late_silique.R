# load lib
library("vcfR")

# import & reformat vcf file 
late_silique_505 <- read.vcfR("~/505/vcf_late_silique_131_sample/combined/505_filtered.vcf.gz")
save(late_silique_505, file="~/505/output/131_sample/late_silique_505.Rdata")

# reform the data 
temp <- late_silique_505

  vcfbi <- is.biallelic(temp) # return with logics indicating whehter biallelic or not...
  vcfref <- subset(getREF(temp), subset = vcfbi) # get ref allele
  vcfalt <- subset(getALT(temp), subset = vcfbi) # get alt allele
  vcfchrom <- subset(getCHROM(temp), subset = vcfbi) # get chromosome info
  vcfpos <- subset(getPOS(temp), subset = vcfbi) # get pos info
  vcfgts <- subset(extract.gt(temp, element = "GT", IDtoRowNames = F), subset = vcfbi)

  temp2 <- data.frame(cbind(vcfchrom,vcfpos,vcfref,vcfalt,vcfgts))
  colnames(temp2)[1:4] <- c("CHROM","POS","REF","ALT")
  rnames <- rownames(temp2)
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/0","0",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/1","1",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("1/1","2",x)))
  row.names(temp2) <- rnames
  rownames(temp2) <- paste(temp2$CHROM, temp2$POS, sep="_")

# GD file
late_silique_505.2 <- temp2[, -(1:4)]
late_silique_505.2.t <- t(late_silique_505.2)
save(late_silique_505.2.t, file = "~/505/output/131_sample/late_silique_505.2.Rdata")
write.table(late_silique_505.2.t, file = "~/505/output/131_sample/late_silique_505.2.txt")

# GM file



# phenotype data
