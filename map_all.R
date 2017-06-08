library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(reshape)
library(dplyr)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.15/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/twopts.f2.05.12.Rdata")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/group.Rdata")

# order marker within each chromosome
mark.f2 <- list()
LG.f2.ord <- list()

for (i in 1:length(group)) {
	mark.f2[[i]] <- make.seq(twopts.f2.05.12, group[[i]])
	LG.f2.ord[[i]] <- order.seq(input.seq = mark.f2[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
}

save(mark.f2, LG.f2.ord, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/LG.f2.ord.Rdata")




