library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(reshape)
library(dplyr)

load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/twopts.f2.04.21.Rdata")

mark.f2.A01 <- make.seq(twopts.f2.04.21, 1:1095)
LG1.f2.ord.A01 <- order.seq(input.seq = mark.f2.A01, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T) 


save(LG1.f2.ord.A01, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/LG1.f2.ord.A01.Rdata")




