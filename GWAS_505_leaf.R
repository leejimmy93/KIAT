# load library
source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("EMMREML") 
install.packages("scatterplot3d")

library(multtest)
library(gplots)
library(LDheatmap) # this could not be succesfully installed 
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

# install GAPIT package 
source("http://zzlab.net/GAPIT/gapit_functions.txt")
# intall EMMA package 
source("http://zzlab.net/GAPIT/emma.txt")

######### load data
# genotype data 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output")
load("~/505/output/leaf_505.Rdata")

# reform the data
temp <- leaf_505 

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
  save(temp2, file="/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/temp2.Rdata") 
  rownames(temp2) <- paste(temp2$CHROM, temp2$POS, sep="_")

# GD file 
leaf_505.2 <- temp2[, 5:98] 
leaf_505.2.t <- t(leaf_505.2)
save(leaf_505.2.t, file = "leaf_505.2.Rdata")
write.table(leaf_505.2.t, file = "leaf_505.2.txt")

# leaf_505.2.matrix <- data.matrix(leaf_505.2.t)
# GM file 
 
SNP.info <- data.frame(matrix(nrow = ncol(leaf_505.2.t), ncol = 3))
colnames(SNP.info) <- c("Name", "Chromosome", "Position")
SNP.info$Name <- colnames(leaf_505.2.t)
SNP.info$Chromosome <- gsub("([[:print:]]+)(_)([[:digit:]]+)", "\\1", colnames(leaf_505.2.t))
SNP.info$Position <- gsub("([[:print:]]+)(_)([[:digit:]]+)", "\\3", colnames(leaf_505.2.t))

# phenotype data 
load("")




















