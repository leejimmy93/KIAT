library(tidyverse)
library(qtl)
library(snowfall) 

sfInit(parallel = TRUE, cpus = parallel::detectCores())

load("~/F2/data/QTL_analysis/LG.f2.after.crossover.Rdata")
LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32)
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1)

sfExport("LG.f2.after.crossover")
sfLibrary(snowfall)

scanone.perm.imp <- 
  ssfExport("LG.f2.after.crossover")
sfLibrary(snowfall)apply(seq_along(LG.f2.after.crossover$pheno), function(trait){
    print(trait)
    tmp <-scanone(LG.f2.after.crossover,pheno.col = trait, method="imp",n.perm=1000)
    summary(tmp)[1]
  })


