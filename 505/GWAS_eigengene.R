library(tidyverse)
load(file = "~/505/network_analysis/output/MEs_with_ID.Rdata")

MEs <- 
MEs %>% 
  mutate(Taxa = ID) %>% 
  dplyr::select(Taxa, MEgreen:MEturquoise) 

library(multtest)
library(gplots)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

setwd("~/505/network_analysis/output/GWAS_eigengene/")

# genotype data
geno_505_hmp <- read.table("~/505/vcf_late_silique_131_sample/combined/505_filtered_het_0.2.recode.hmp.txt", head=FALSE)
geno_505_hmp$V3 <- as.numeric(as.factor(geno_505_hmp$V3))
geno_505_hmp[1,3] <- "chrom"

colnames(MEs)[40] <- "Taxa"
pheno_data_505 <- MEs

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
