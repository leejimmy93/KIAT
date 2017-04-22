# load package 
# install.packages("onemap")
# install.packages("tkrplot", type="source")
library(tcltk)
library(tkrplot)
library(onemap)
library(ggplot2)
library(reshape)

# load data 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output")
F2.data <- read.mapmaker(file="F2_geno_for_one_map_final.txt")

# estimate two-point rf
twopts.f2.04.21 <- rf.2pts(F2.data, LOD = 6, max.rf = 0.25)

# save & export data 
save(twopts.f2.04.21, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/twopts.f2.04.21Rdata')
