# GWAS result plot 
library("qqman")

Oil_content <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Oil_content.GWAS.Results.csv", header = T)

# Caprylic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Caprylic_acid.GWAS.Results.csv", header = T)

# Capric_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Capric_acid.GWAS.Results.csv", header = T)

# Lauric_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Lauric_acid.GWAS.Results.csv", header = T)

Myristic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Myristic_acid.GWAS.Results.csv", header = T)

Pentadecanoic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Pentadecanoic_acid.GWAS.Results.csv", header = T)

Palmitic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Palmitic_acid.GWAS.Results.csv", header = T)

Palmitoliec_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Palmitoliec_aicd.GWAS.Results.csv", header = T)

Heptadecanoic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Heptadecanoic_acid.GWAS.Results.csv", header = T)

Stearic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Stearic_acid.GWAS.Results.csv", header = T)

Oleic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Oleic_acid.GWAS.Results.csv", header = T)

vaccenic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..vaccenic_acid.GWAS.Results.csv", header = T)

Linoleic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Linoleic_acid.GWAS.Results.csv", header = T)

Arachidic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Arachidic_acid.GWAS.Results.csv", header = T)

cis_11_Eicosenoic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..cis_11_Eicosenoic_acid.GWAS.Results.csv", header = T)

Linolenic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Linolenic_acid.GWAS.Results.csv", header = T)

Behenic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Behenic_acid.GWAS.Results.csv", header = T)

Erucic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/fatty_acid_all/GAPIT..Erucic_acid.GWAS.Results.csv", header = T)

colnames(Oil_content)[1:4] <- c("SNP", "CHR", "BP", "P")
colnames(Oil_content)[1:4] <- c("SNP", "CHR", "BP", "P")
colnames(bolting_time)[1:4] <- c("SNP", "CHR", "BP", "P")


png("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/figure/131_sample/bolting_flowering_oil.png", width=12, height=10, units="in", res=300)

par(mfrow=c(3,1)) 
manhattan(flower_time, main = "Flowering time", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), supggestiveline = 5, genomewideline = F 
    ) 

manhattan(Oil_content, main = "Oil content", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = 5, genomewideline = F 
    ) 

manhattan(bolting_time, main = "bolting_time", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = 5, genomewideline = F
    )
 
dev.off() 
