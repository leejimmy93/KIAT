# load library
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("multtest")
# install.packages("gplots")
# install.packages("LDheatmap")
# install.packages("genetics")
# install.packages("EMMREML") 
# install.packages("scatterplot3d")

library(multtest)
library(gplots)
# library(LDheatmap) # this could not be succesfully installed 
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

# install GAPIT package 
source("http://zzlab.net/GAPIT/gapit_functions.txt")
# intall EMMA package 
source("http://zzlab.net/GAPIT/emma.txt")

# set working directory 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/GWAS/tmp2")

######### load data
# genotype data 
geno_505_hmp <- read.table("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/WGS/vcf_raw/filtered/filter_1/hmp/505.hmp.txt", head=FALSE, stringsAsFactors = F)
# get chrom name to numerics
geno_505_hmp[1,3] <- "chrom"
geno_505_hmp[1,] <- gsub("_S.*", "", unlist(geno_505_hmp[1,])) 
geno_505_hmp$V3 <- as.numeric(as.factor(geno_505_hmp$V3))

# phenotype data 
load("~/505/WGS/GWAS/fatty_acid_oil_K_all_sub.Rdata")
pheno_data_505 <- fatty_acid_oil_K_all_sub

myY <- pheno_data_505
myG <- geno_505_hmp

# run GAPIT 
myGAPIT <- GAPIT( 
Y=myY,
G=myG, 
PCA.total=4,
Geno.View.output=FALSE, 
PCA.View.output=FALSE,
Model.selection = TRUE 
)  

print("finished, YAAH!")


