library(qtl)
# library(snowfall)

load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/QTL_analysis/LG.f2.after.crossover_fatty_acid_lipid.Rdata")

LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32) # imputation? 
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1) 

# sfInit(parallel = TRUE, cpus = 16) 
# sfExport("LG.f2.after.crossover")
# sfLibrary(qtl)

# system.time(
# scanone.perm.imp <- 
#   sfLapply(seq_along(LG.f2.after.crossover$pheno), function(trait){
#     print(trait) # print doesn't work in here 
#     tmp <-scanone(LG.f2.after.crossover,pheno.col = trait, method="imp",n.perm=1000)
#     summary(tmp)[1] # #keep the 95th percentile for future use.This corresponds to p <0.05
#   }) # takes 40 mins to finish 
# )
# sfStop() 

# names(scanone.perm.imp) <- colnames(LG.f2.after.crossover$pheno) 

system.time(
scanone.imp <- 
lapply(seq_along(LG.f2.after.crossover$pheno), function(trait) {
  print(trait)
  scanone(LG.f2.after.crossover,pheno.col=trait,method="imp")
}) 
)
names(scanone.imp) <- colnames(LG.f2.after.crossover$pheno) 

save(scanone.imp, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/QTL_analysis/scanone.imp_fatty_acid_lipid.Rdata")
