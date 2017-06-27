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
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample_SN")

######### load data
# genotype data 
geno_505_hmp <- read.table("~/505/leaf_VS_late_silique/geno_505_131sample.hmp.txt", head=FALSE)
# get chrom name to numerics
geno_505_hmp$V3 <- as.numeric(geno_505_hmp$V3)
geno_505_hmp[1,3] <- "chrom"

# phenotype data 
load("~/505/data/pheno_data_505.silique.SN.Rdata")
pheno_data_505 <- pheno_data_505.silique.SN 

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

print("finished, YAAH!")


