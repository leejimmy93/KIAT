library(qtl)

LG.f2 <- read.cross("mm", file = "~/F2/data/QTL_analysis/F2_geno_for_one_map_final_all_expressed_gene_no_scale_center.txt", mapfile = "~/F2/data/QTL_analysis/LG.f2.madmapper.map") ### takes one day to read in this data...

save(LG.f2, file = "~/F2/output/QTL_analysis/LG.f2.all_expressed_genes_no_scale_center.Rdata") 
