library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(reshape)
library(dplyr)

F2.data <- read.mapmaker(file="~/F2/output/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/twopts.f2.04.21.Rdata")

mark.f2.A09 <- make.seq(twopts.f2.04.24, 2839:3392)
LG1.f2.ord.A09 <- order.seq(input.seq = mark.f2.A09, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T) 


save(LG1.f2.ord.A09, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/LG1.f2.ord.A09.Rdata")




