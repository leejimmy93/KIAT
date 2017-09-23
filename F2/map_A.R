library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
# library(reshape)
library(dplyr)
library(snowfall)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.10/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/twopts.LOD3_rf0.5.Rdata")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/group.Rdata")

# order marker within each chromosome
sfInit(parallel = TRUE, cpus = 12)
sfLibrary(onemap)
sfExport("F2.data", "twopts.f2.LOD3_rf0.5", "group.A")

mark.f2.A <- list()
LG.f2.ord.A <- list()

for (i in 1:length(group.A)) {
	mark.f2.A[[i]] <- make.seq(twopts.f2.LOD3_rf0.5, group.A[[i]])
	LG.f2.ord.A[[i]] <- order.seq(input.seq = mark.f2.A[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
}

sfStop()

save(mark.f2.A, LG.f2.ord.A, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/LG.f2.ord.A.Rdata")




