# use 0.55 as max.rf to capture more markers in case...  
# load package 
# install.packages("onemap")
# install.packages("tkrplot", type="source")
library(tcltk)
library(tkrplot)
library(onemap)
library(ggplot2)

# load data 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10")
F2.data <- read.mapmaker(file="F2_geno_for_one_map_final.txt")

# estimate two-point rf
twopts.f2.LOD4_rf0.5 <- rf.2pts(F2.data, LOD = 4, max.rf = 0.5)

# save & export data 
save(twopts.f2.LOD4_rf0.5, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/twopts.LOD4_rf0.5.Rdata')





