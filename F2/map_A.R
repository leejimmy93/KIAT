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
mark.f2.A <- list()
LG.f2.ord.A <- list()

for (i in 1:length(group.A)) {
	mark.f2.A[[i]] <- make.seq(twopts.f2.05.12, group.A[[i]])
	LG.f2.ord.A[[i]] <- order.seq(input.seq = mark.f2.A[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
}

save(mark.f2.A, LG.f2.ord.A, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.15/LG.f2.ord.A.Rdata")




