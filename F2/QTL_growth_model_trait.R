library(qtl)
library(tidyverse)
library(snowfall)

cross.F2 <- read.cross("csvsr", genfile =  "~/F2/data/QTL_analysis/LG.f2.madmapper.final_gen_revised.csv", 
                         phefile = "~/F2/output/pheno/growth_model_trait.2.csv",
                         genotypes = c("AA", "AB", "BB"))   # although the sample IDs are not matched in the original phe and gen file, I still get the right result. 

# cross.F2$pheno %>% colnames() # for some reason, id is always read as a phenotype 

cross.F2 <- sim.geno(cross.F2,step=1,n.draws=32) # imputation?  
cross.F2 <- calc.genoprob(cross.F2,step=1) 

cim_growth_model_trait.F2 <- 
lapply(seq_along(cross.F2$pheno[1:12]), function(trait) {
  print(trait)
  cim(cross.F2, n.marcovar=5, pheno.col=trait,method="em")
})  

names(cim_growth_model_trait.F2) <- colnames(cross.F2$pheno[1:12])

# scanone_growth_model_trait.F2 %>% colnames()
# plot(scanone_growth_model.F2[,c(1,2,3)])

sfInit(parallel = TRUE, cpus = 4)  
sfExport("cross.F2")
sfLibrary(qtl)
 
cim.perm <- 
  sfLapply(seq_along(cross.F2$pheno[1:12]), function(trait){
    message(trait) # message doesn't work in here either 
    tmp <- cim(cross.F2,
               pheno.col = trait, 
               n.marcovar=5, 
               method="em",
               n.perm=1000)
    summary(tmp)[1] # #keep the 95th percentile for future use.This corresponds to p <0.05
  }) # takes almost 4 hours to finish 

sfStop() 

names(cim.perm) <- colnames(cross.F2$pheno[1:12])

save(cim_growth_model_trait.F2, cim.perm, file = "~/F2/output/growth_model/cim_growth_model_trait.Rdata")



