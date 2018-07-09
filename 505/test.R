library(tidyverse)
load(file = "~/505/network_analysis/output/MEs_with_ID.Rdata")
load("~/505/network_analysis/output/pheno_data_505.Rdata")
# genotype data
geno_505_hmp <- read.table("~/505/vcf_late_silique_131_sample/combined/505_filtered_het_0.2.recode.hmp.txt", head=FALSE)
geno_505_hmp$V3 <- as.numeric(geno_505_hmp$V3) # don't use as.factor() here
geno_505_hmp[1,3] <- "chrom"

MEs <- 
MEs %>% 
  mutate(Taxa = ID) %>% 
  dplyr::select(Taxa, MEgreen:MEturquoise) 

ID.2 <- gsub("\\_|\\-", "", MEs$Taxa)
ID.2 <- gsub("(505)(K)([[:digit:]]+)", "ID_\\1_\\2\\3", ID.2)

ID.matching <- data.frame(ID = MEs$Taxa,
                          ID.2 = ID.2)

pheno_data_505 <- 
pheno_data_505 %>% 
  left_join(ID.matching, by=c("ID" = "ID.2")) %>% 
  mutate(Taxa = ID.y) %>%
  left_join(MEs) %>% 
  dplyr::select(Taxa, Erucic_acid_year2, MEsalmon) 

# plot(pheno_data_505$Erucic_acid_year2, pheno_data_505$MEsalmon)

library(multtest) 
library(gplots)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

myY <- pheno_data_505 
myG <- geno_505_hmp 

setwd("~/505/network_analysis/output/test/") 
# run GAPIT
myGAPIT <- GAPIT(
Y=myY,
G=myG,
PCA.total=0,
Geno.View.output=FALSE,
PCA.View.output=FALSE,
Model.selection = TRUE
)  
