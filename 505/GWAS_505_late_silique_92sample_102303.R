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
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_92_102303")

######### load data
# genotype data 
geno_505_hmp <- read.table("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/leaf_VS_late_silique/505_filtered_common_late_silique_92.hmp.txt", head=FALSE)
# get chrom name to numerics
geno_505_hmp$V3 <- as.numeric(geno_505_hmp$V3)
geno_505_hmp[1,3] <- "chrom"

# phenotype data 
load("~/505/data/phenotype/gh.seed.bolt.table.wide.Rdata")
# reform phentype data to have match taxa ID with geno data
# get taxa ID from sample description file 
sample_des_c <- read.csv("~/505/data/phenotype/batch_c.csv", header=T)
sample_des_d <- read.csv("~/505/data/phenotype/batch_d.csv", header=T)
sample_des_e <- read.csv("~/505/data/phenotype/batch_e.csv", header=T)

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

# get the 92 common lines between leaf & late silique
load("~/505/output/common.lines.leaf.131.sample.Rdata")
pheno_data_505_92 <- pheno_data_505[which(pheno_data_505$Taxa %in% common.lines$Taxa.y),]
rownames(pheno_data_505_92) <- 1:nrow(pheno_data_505_92)

myY <- pheno_data_505_92
myG <- geno_505_hmp

# run GAPIT 
myGAPIT <- GAPIT( 
Y=myY,
G=myG, 
PCA.total=3,
Geno.View.output=FALSE, 
PCA.View.output=FALSE 
)  

print("finished, YAAH!")


