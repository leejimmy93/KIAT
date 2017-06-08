# The reason to try default LOD3 & max.rf=0.5 is to see whether I can get differnt linkage groups use these parameters. 
# load package 
# install.packages("onemap")
# install.packages("tkrplot", type="source")
library(tcltk)
library(tkrplot)
library(onemap)
library(ggplot2)
library(reshape)

# load data 
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15")
F2.data <- read.mapmaker(file="F2_geno_for_one_map_final.txt")

# estimate two-point rf
twopts.f2.LOD3_rf0.5 <- rf.2pts(F2.data)

# save & export data 
save(twopts.f2.LOD3_rf0.5, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/twopts.LOD3_rf0.5.Rdata')

# assign marker to linkage group # this way we get all the markers to be assigned to one linkage group, bad...
mark.all.f2.LOD3_rf0.5 <- make.seq(twopts.f2.LOD3_rf0.5, "all")
LGs.f2.LOD3_rf0.5 <- group(mark.all.f2.LOD3_rf0.5)
save(LGs.f2.LOD3_rf0.5, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/LGs.f2.LOD3_rf0.5.Rdata')

# Julin suggested to assign linkage groups based on their mapping position on the chromosome
# get the number of markers on each chromosome



