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
twopts.f2.04.24 <- rf.2pts(F2.data, LOD = 6, max.rf = 0.25)

# save & export data 
save(twopts.f2.04.24, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/twopts.f2.04.21.Rdata')

# assign marker to linkage group # this way we get all the markers to be assigned to one linkage group, bad...
# mark.all.f2 <- make.seq(twopts.f2.04.21, "all")
# LGs.f2 <- group(mark.all.f2, LOD = 6, max.rf = 0.25)
# save(LGs.f2, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/LGs.f2.04.24.Rdata')

# Julin suggested to assign linkage groups based on their mapping position on the chromosome
# get the number of markers on each chromosome



