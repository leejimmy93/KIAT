library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(reshape)
library(dplyr)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.15/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/twopts.f2.LOD3_rf0.5.Rdata")

mark.f2.C06 <- make.seq(twopts.f2.LOD3_rf0.5, 2644:2728)
LG.f2.ord.C06 <- order.seq(input.seq = mark.f2.C06, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T) 


save(LG.f2.ord.C06, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/LG.f2.ord.C06.Rdata")




