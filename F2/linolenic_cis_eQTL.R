library(tidyverse)
library(qtl)

load("~/F2/output/eQTL/cis_trans_result_new_flipped.Rdata") 

allele.effect <- function(gene_ID){
  cross.F2 <- read.cross("csvsr", 
                         genfile ="~/F2/data/QTL_analysis/LG.f2.madmapper.final_gen_revised_flipped.csv", 
                         phefile = "~/F2/output/network_analysis/vstMat.f2.batch.corrected_revised.csv",
                         genotypes = c("AA", "AB", "BB"))   # although the sample IDs are not matched in the original phe and gen file, I still get the right result. 
  cross.F2$pheno <- as.data.frame(cross.F2$pheno[,gene_ID]) 
  cross.F2 <- sim.geno(cross.F2,step=1,n.draws=32) # imputation?  
  cross.F2 <- calc.genoprob(cross.F2,step=1) 

  chr = cis_eQTL[cis_eQTL$gene_ID == gene_ID,]$eQTL_chr
  chr = gsub("chr", "", chr)
  pos = cis_eQTL[cis_eQTL$gene_ID == gene_ID,]$pos
  qtlm <- makeqtl(cross.F2,chr=chr,pos=pos) # Ruijuan: QTL with interactions

#make a model for QTL action.  Since there is no evidence of interaction, start with an additive model
  qtl.fit <- fitqtl(cross=cross.F2,qtl=qtlm,formula = y ~ Q1, get.ests=T)

#examine our model
  summary(qtl.fit)    
} 

cis_eQTL.qtl.combined.final <- read.csv("~/F2/output/eQTL/cis_eQTL.qtl.combined.final.csv", row.names = 1) 

cis_eQTL.qtl.combined.final.linolenic <- 
cis_eQTL.qtl.combined.final %>% 
  dplyr::filter(trait.y == "Linolenic_acid") 

system.time(
cis_eQTL_effect <-   
lapply(cis_eQTL.qtl.combined.final.linolenic$gene_ID %>% as.character(), function(gene) {
  allele.effect(gene)
}
))

save(cis_eQTL_effect, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/for_paper/linolenic_acid/cis_eQTL_effect.Rdata")


