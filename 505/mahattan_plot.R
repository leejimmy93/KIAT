# GWAS result plot 
library("qqman")

Erucic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/filter_by_het/GAPIT..Erucic_acid.GWAS.Results.csv", header = T)

Oil_content <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/filter_by_het/GAPIT..Oil_content.GWAS.Results.csv", header = T)

Oleic_acid <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/filter_by_het/GAPIT..Oleic_acid.GWAS.Results.csv", header = T)

colnames(Erucic_acid)[1:4] <- c("SNP", "CHR", "BP", "P")
colnames(Oil_content)[1:4] <- c("SNP", "CHR", "BP", "P")
colnames(Oleic_acid)[1:4] <- c("SNP", "CHR", "BP", "P")

png("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/figure/131_sample/fatty_acid.png", width=12, height=10, units="in", res=300)

par(mfrow=c(3,1)) 
manhattan(Erucic_acid, main = "Erucic acid", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), supggestiveline = 5, genomewideline = F 
    ) 

manhattan(Oleic_acid, main = "Oleic acid", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = 5, genomewideline = F 
    ) 

manhattan(Oil_content, main = "Oil content", ylim = c(0, 8), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = 5, genomewideline = F
    )
 
dev.off() 
