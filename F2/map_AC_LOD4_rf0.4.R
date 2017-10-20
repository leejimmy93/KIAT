library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(dplyr)
library(snowfall)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.10/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/twopts.LOD4_rf0.4.Rdata")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/group.AC.Rdata")

mark.f2 <- list()
LG.f2.ord <- list()

system.time(
for (i in 1:length(group.AC)) {
	mark.f2[[i]] <- make.seq(twopts.f2.LOD4_rf0.4, group.AC[[i]])
	LG.f2.ord[[i]] <- order.seq(input.seq = mark.f2[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 4,
                        draw.try = T, wait = 1
                        )
}
)

save(mark.f2, LG.f2.ord, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/LG.f2.ord.AC.LOD4_rf0.4.Rdata")




