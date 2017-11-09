library(tidyverse)
library(qtl)
library(snowfall) 

# already finished 1000 permutation 

# run QTL anlaysis on all 43 traits (scanone)
load("~/F2/output/QTL_analysis/LG.f2.after.crossover.43traits.Rdata")
summaryMap(LG.f2.after.crossover)
LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32) # imputation? 
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1) # calculate the probability of the true underlying genotypes given the observed multipoint marker data --> for each imputed data, give a probability?  

print("scanone QTL")

system.time(
scanone.imp <- 
lapply(seq_along(LG.f2.after.crossover$pheno), function(trait) {
  print(trait)
  scanone(LG.f2.after.crossover,pheno.col=trait,method="imp")
}) 
)
names(scanone.imp) <- colnames(LG.f2.after.crossover$pheno)

save(scanone.imp, file = "~/F2/output/QTL_analysis/scanone.imp.43traits")

# run cim mapping permutation and save the result 
print("cim QTL") 

cim.qtl <- 
lapply(seq_along(LG.f2.after.crossover$pheno), function(trait) {
  print(trait)
  cim(LG.f2.after.crossover, n.marcovar=5, pheno.col=trait,method="em")
}) 

names(cim.qtl) <- colnames(LG.f2.after.crossover$pheno)

print("cim permutation")

sfInit(parallel = TRUE, cpus = parallel::detectCores()) 
sfExport("LG.f2.after.crossover")
sfLibrary(qtl)  

# permutation 
cim.perm <- 
  sfLapply(seq_along(LG.f2.after.crossover$pheno), function(trait){
    message(trait) # message doesn't work in here either 
    tmp <- cim(LG.f2.after.crossover,
               pheno.col = trait, 
               n.marcovar=5, 
               method="em",
               n.perm=1000)
    summary(tmp)[1] # #keep the 95th percentile for future use.This corresponds to p <0.05
  }) # takes almost 4 hours to finish 

sfStop() 

names(cim.perm) <- colnames(LG.f2.after.crossover$pheno)

# save output 
save(cim.qtl, file = "~/F2/output/QTL_analysis/cim.qtl.43traits.Rdata")
save(cim.perm, file = "~/F2/output/QTL_analysis/cim.perm.43traits.Rdata")
