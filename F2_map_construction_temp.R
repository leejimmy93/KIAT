# load package 
install.packages("onemap")
install.packages("tkrplot", type="source")
library(tcltk)
library(tkrplot)
library(onemap)
library(ggplot2)
library(reshape)

# load data 
load("~/F2/data/data.Rdata")
# estimate two-point rf
twopts.f2 <- rf.2pts(data, LOD = 6, max.rf = 0.25)
# save & export data 
save(twopts.f2, file='/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/twopts.f2.Rdata')
