library(qtl)

load("~/F2/output/QTL_analysis/cim.qtl.43traits.Rdata")
load("~/F2/output/QTL_analysis/LG.f2.after.crossover.43traits.Rdata")

cim.perm.3.traits <- 
  lapply(c(5, 12, 14), function(trait){
    print(trait) # message doesn't work in here either 
    tmp <- cim(LG.f2.after.crossover, 
               pheno.col = trait, 
               n.marcovar=5, 
               method="em",
               n.perm=1000)
    summary(tmp)[1] # #keep the 95th percentile for future use.This corresponds to p <0.05
  }) # takes almost 4 hours to finish 

names(cim.perm.3.traits) <- colnames(LG.f2.after.crossover$pheno)

save(cim.perm.3.traits, file="~/F2/output/QTL_analysis/cim.perm.3.traits.Rdata") 
