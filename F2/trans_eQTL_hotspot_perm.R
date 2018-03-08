library(tidyverse)
library(qtl)
library(qtlhot)
library(snowfall) 

load("~/F2/output/eQTL/scanone-eqtl_F2_flipped.RData")
scanone_eQTL.F2 %>% dim() # 4887 56182 

scanone_eQTL.F2[1:10, 1:10]

### check correlation, later... 
cross.F2
lod.thrs
alphas <- seq(0.01, 0.10, by=0.01)
alphas

### LOD threshold level for e-trait
lod.thr <- lod.thrs[5]
lod.thr 

### get only significant intervals for each e-trait, using LOD drop method 
high1 <- highlod(scanone_eQTL.F2, lod.thr = min(lod.thrs), drop.lod = 1.5)
max(high1, lod.thr = lod.thrs) # max number of e-trait fall into loci with different lod threshold 

hots1 <- hotsize(high1, lod.thr = lod.thr) 
summary(hots1) # for each genomic position, 

# plot(hots1, cex.lab = 1.5, cex.axis = 1.5) 

### permutation to get the statistical significance # permutation takes time, tried to 
set.seed(12345) 

cross.F2$pheno <- cross.F2$pheno[,-1]

system.time(
hotperm1 <- hotperm(cross = cross.F2, 
n.quant = 600,
n.perm = 100,
lod.thrs = lod.thrs,
alpha.levels = alphas,
drop.lod = 1.5,
verbose = FALSE)
)   

save(hotperm1, file = "~/F2/output/eQTL/hotperm1.Rdata")
