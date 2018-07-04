library(multtest)
library(gplots)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
library(snowfall)

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

load("~/505/network_analysis/output/GWAS_trait_input.Rdata")
setwd("~/505/network_analysis/output/GWAS_traits/")

myY <- myY[[4]]
# run GAPIT
myGAPIT <- GAPIT(
    Y=myY,
    G=myG,
    PCA.total=0,
    Geno.View.output=FALSE,
    PCA.View.output=FALSE,
    Model.selection = TRUE
    ) 

