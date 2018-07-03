load("~/505/network_analysis/output/pheno_data_505.Rdata")

library(multtest)
library(gplots)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
library(tidyverse)

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

setwd("~/505/network_analysis/output/GWAS_traits/")

# genotype data
geno_505_hmp <- read.table("~/505/vcf_late_silique_131_sample/combined/505_filtered_het_0.2.recode.hmp.txt", head=FALSE)
geno_505_hmp$V3 <- as.numeric(as.factor(geno_505_hmp$V3))
geno_505_hmp[1,3] <- "chrom"

geno_505_hmp[1,] <- gsub("\\_|\\-", "", geno_505_hmp[1,])
geno_505_hmp[1,] <- gsub("(505)(K)([[:digit:]]+)", "ID_\\1_\\2\\3", geno_505_hmp[1,])

colnames(pheno_data_505)[1] <- "Taxa"

myY <- pheno_data_505 
myG <- geno_505_hmp 

# run GAPIT
myGAPIT <- GAPIT(
Y=myY,
G=myG,
PCA.total=0,
Geno.View.output=FALSE,
PCA.View.output=FALSE,
Model.selection = TRUE
)    
