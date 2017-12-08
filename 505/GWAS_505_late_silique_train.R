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
library(tidyverse)

# install GAPIT package 
source("http://zzlab.net/GAPIT/gapit_functions.txt")
# intall EMMA package 
source("http://zzlab.net/GAPIT/emma.txt")

# set working directory 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/genomic_prediction/output")

######### load data
# genotype data 
geno_505_hmp <- read.table("~/505/vcf_late_silique_131_sample/combined/505_filtered_het_0.2.recode.hmp.txt", head=FALSE, stringsAsFactors=FALSE)
# get chrom name to numerics
geno_505_hmp$V3 <- as.numeric(as.factor(geno_505_hmp$V3))
geno_505_hmp[1,3] <- "chrom"
# geno_505_hmp <- geno_505_hmp[,-91] # remove K88 because it is empty, need to work on it later
# subset the 85 training samples
load("~/505/genomic_prediction/output/train_sample.Rdata")
more <- geno_505_hmp[1,1:11] %>% as.character()
to_keep <- c(more, train_sample)
geno_505_hmp <- geno_505_hmp[,geno_505_hmp[1,] %in% to_keep]
geno_505_hmp[]

# phenotype data 
load("~/505/data/phenotype/gh.seed.bolt.table.wide.Rdata")
# reform phentype data to have match taxa ID with geno data
# get taxa ID from sample description file 
sample_des_c <- read.csv("~/505/data/phenotype/batch_c.csv", stringsAsFactors=F, header=T)
sample_des_d <- read.csv("~/505/data/phenotype/batch_d.csv", stringsAsFactors=F, header=T)
sample_des_e <- read.csv("~/505/data/phenotype/batch_e.csv", stringsAsFactors=F, header=T)

sample_des_c_sub2 <- sample_des_c[,c("Sample.ID", "No..of.Sowing", "Name")]
sample_des_d_sub2 <- sample_des_d[,c("Sample.ID", "No..of.Sowing", "Name")]
sample_des_e_sub2 <- sample_des_e[,c("Sample.ID", "No..of.Sowing", "Name")]

sample_des_d_sub2 <- sample_des_d_sub2[c(1:41),]
sample_des_sub2 <- rbind(sample_des_c_sub2, sample_des_d_sub2, sample_des_e_sub2)

# add name to to replace SN 
sample_des_sub2$SN <- paste("K", sample_des_sub2$No..of.Sowing, sep="")
gh.seed.bolt.table.wide.merge <- merge(sample_des_sub2, gh.seed.bolt.table.wide, by = "SN")
pheno_data_505 <- gh.seed.bolt.table.wide.merge[,c("Sample.ID", "FATTY ACID_Erucic acid (C22:1n9)", "FATTY ACID_oil contents (%)", "FATTY ACID_Oleic acid (C18:1n9c)")]
colnames(pheno_data_505) <- c("Taxa", "Erucic_acid", "Oil_content", "Oleic_acid")
pheno_data_505
K250_K251 <- data.frame(Taxa = c("505_K_250", "505_K_251"),
			Erucic_acid = rep(NA, 2),
			Oil_content = rep(NA, 2),
			Oleic_acid = rep(NA, 2))
pheno_data_505 <- rbind(pheno_data_505, K250_K251)
 
# subset the 85 training samples 
pheno_data_505 <- pheno_data_505[pheno_data_505$Taxa %in% train_sample,]
# pheno_data_505 <- pheno_data_505[-81,] 

myY <- pheno_data_505
myG <- geno_505_hmp

# subset to get only 1st chromosome result 
# myG.01 <- myG[myG$V3==1,]
# myG.01.final <- data.frame(rbind(myG[1,], myG.01)) 

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


