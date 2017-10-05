library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(snowfall)
library(dplyr)

F2.data <- read.mapmaker(file="~/F2/output/missing_rate_0.10/F2_geno_for_one_map_final.txt")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/twopts.LOD3_rf0.55.Rdata")
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/group.AC.Rdata")

# order marker within each chromosome
sfInit(parallel = TRUE, cpus = 12)
sfLibrary(onemap)
sfExport("F2.data", "twopts.f2.LOD3_rf0.55", "group.AC")

LG.f2.ord.0.55 <-
sfLapply(1:length(group.AC), function(LG) {
        tmp <- make.seq(twopts.f2.LOD3_rf0.55, group.AC[[LG]])
        	order.seq(input.seq = tmp, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
})

sfStop()

save(LG.f2.ord.0.55, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/missing_rate_0.10/LG.f2.ord.AC.0.55.Rdata")




