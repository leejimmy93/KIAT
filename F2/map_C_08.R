library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
# library(reshape)
library(dplyr)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.10/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/twopts.LOD3_rf0.5.Rdata")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/group.Rdata")

# order marker within each chromosome
mark.f2.C8 <- list()
LG.f2.ord.C8 <- list()

for (i in 10) {
	mark.f2.C8[[i]] <- make.seq(twopts.f2.LOD3_rf0.5, group.C[[i]])
	LG.f2.ord.C8[[i]] <- order.seq(input.seq = mark.f2.C8[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
}

save(mark.f2.C8, LG.f2.ord.C8, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/LG.f2.ord.C8.Rdata")




