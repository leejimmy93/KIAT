library(qtl)

load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/QTL_analysis/LG.f2.after.crossover_Erucic_Oleic.Rdata")

LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32)
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1)

system.time(
scantwo.imp <-
lapply(seq_along(LG.f2.after.crossover$pheno), function(trait) {
  print(trait)
  scantwo(LG.f2.after.crossover,pheno.col=trait,method="imp")
})
)
names(scantwo.imp) <- colnames(LG.f2.after.crossover$pheno)

# save output
save(scantwo.imp, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/scantwo/scantwo.imp.Erucic_Oleic.Rdata")
