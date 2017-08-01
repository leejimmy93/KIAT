# GWAS result plot 
library("qqman")

flower_time <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/bolting_floweirng_fruition_oil/GAPIT..flowering_time.GWAS.Results.csv", header = T)

Oil_content <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/bolting_floweirng_fruition_oil/GAPIT..oil_content.GWAS.Results.csv", header = T)

bolting_time <- read.csv("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/505/output/myGAPIT/late_silique_131_sample/bolting_floweirng_fruition_oil/GAPIT..bolting_time.GWAS.Results.csv", header = T)

colnames(flower_time)[1:4] <- c("SNP", "CHR", "BP", "P")
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
