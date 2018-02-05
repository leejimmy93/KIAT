library(tidyverse)
library(qtl)

LG.f2.after.crossover <- read.cross("csvsr", genfile = "~/F2/data/QTL_analysis/LG.f2.madmapper.final_gen_revised_flipped.csv",
                     phefile = "~/F2/data/QTL_analysis/F2.pheno.csv",
                     genotypes = c("AA", "AB", "BB")) # the only problem is that 44 phenotypes were read instead of 43, need to figure out why later

LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32) # imputation
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1) # calculate the probability of the true underlying genotypes given the observed multipoint marker data --> for each imputed data, give a probability?

# scanone_Myristic <-
# scanone(LG.f2.after.crossover, pheno.col = "Myristic_acid", method = "em", model = "2part", upper = F)

perm_Myristic <- 
scanone(LG.f2.after.crossover, pheno.col = "Myristic_acid", method = "em", model = "2part", upper = F, n.perm = 1000)

LG.f2.after.crossover$pheno$Heptadecanoic_acid_new <- ifelse(LG.f2.after.crossover$pheno$Heptadecanoic_acid == 0, 0, 1) 

# scanone_Heptadecanoic <-
# scanone(LG.f2.after.crossover, pheno.col = "Heptadecanoic_acid_new", method = "em", model = "binary") # binary
 
perm_Heptadecanoic <-
scanone(LG.f2.after.crossover, pheno.col = "Heptadecanoic_acid_new", method = "em", model = "binary", n.perm = 1000)

save(perm_Myristic, perm_Heptadecanoic, file="~/F2/output/QTL_analysis/perm_new_model.Rdata")
