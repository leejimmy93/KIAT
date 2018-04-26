library(tidyverse)
library(qtl)
library(snowfall)
 
cis_eQTL.qtl.combined.final <- read.csv("~/F2/output/eQTL/cis_eQTL.qtl.combined.final.csv", row.names = 1)

sfInit(parallel = TRUE, cpus = 4)
sfExport("cis_eQTL.qtl.combined.final")
sfLibrary(qtl)
sfSource("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/KIAT/F2/function_cis_eQTL_effect.R")  

system.time(
cis_eQTL_effect_all <-   
sfLapply(unique(cis_eQTL.qtl.combined.final$gene_ID) %>% as.character(), function(gene) {
  allele.effect(gene)
}
))

sfStop()
save(cis_eQTL_effect_all, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/cis_eQTL_effect_all.Rdata")


