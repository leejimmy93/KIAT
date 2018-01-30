library(qtl)
library(tidyverse)

cross.F2 <- read.cross("csvsr", genfile =  "~/F2/data/QTL_analysis/LG.f2.madmapper.final_gen_revised.csv", 
                         phefile = "~/F2/output/pheno/growth_model_trait.csv",
                         genotypes = c("AA", "AB", "BB"))   # although the sample IDs are not matched in the original phe and gen file, I still get the right result. 

# cross.F2$pheno %>% colnames() # for some reason, id is always read as a phenotype 

cross.F2 <- sim.geno(cross.F2,step=1,n.draws=32) # imputation?  
cross.F2 <- calc.genoprob(cross.F2,step=1) 

scanone_growth_model_trait.F2 <- scanone(cross.F2, pheno.col = 2:ncol(cross.F2$pheno), method = "imp", use = "all.obs") 

# scanone_growth_model_trait.F2 %>% colnames()
# plot(scanone_growth_model.F2[,c(1,2,3)])

set.seed(12345)  
system.time(
  permtest.F2 <- scanone(cross.F2, pheno.col = 2:ncol(cross.F2$pheno), method = "imp", n.perm = 1000)  
)

alphas <- seq(0.01, 0.10, by = 0.01) 
lod.thrs <- summary(permtest.F2, alphas) 
lod.thrs

save(scanone_growth_model_trait.F2, lod.thrs, "~/F2/output/growth_model/scanone_growth_model_trait.Rdata")

cim_growth_model_trait.F2 <- cim(cross.F2, n.marcovar=5, pheno.col = 2:ncol(cross.F2$pheno), method = "em")  

set.seed(12345)
system.time(
  permtest.F2.cim <- cim(cross.F2, pheno.col = 2:ncol(cross.F2$pheno), n.marcovar=5, method = "em", n.perm = 1000)
)

alphas <- seq(0.01, 0.10, by = 0.01)
lod.thrs.cim <- summary(permtest.F2.cim, alphas)
lod.thrs.cim

save(cim_growth_model_trait.F2, lod.thrs.cim, "~/F2/output/growth_model/cim_growth_model_trait.Rdata")



