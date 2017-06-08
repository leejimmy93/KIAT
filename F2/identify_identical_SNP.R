# the goal of this script is to remove identical SNPs based on correlation, if pairwise correlation is 1, extract the name of identical SNPs, output file are the SNPs that should be excluded from the dataset (rem.Rdata) and duplicated pairs (dup.pair.Rdata)  
setwd("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output")
load("F2_SNP_correlation.Rdata")

# find duplicate coordinate
dup.cordinate <- which(F2_SNP_correlation >=0.9 & lower.tri(F2_SNP_correlation), arr.ind = T, useNames = F)
dup.cordinate.df <- as.data.frame(dup.cordinate)
sample.ID <- colnames(F2_SNP_correlation)

# extract duplicate pair information based on their coordicate
dup.pair <- data.frame(matrix(nrow = nrow(dup.cordinate), ncol = 2))
for (i in 1:nrow(dup.cordinate)){
 dup.pair[i,1] <- sample.ID[dup.cordinate[i,1]]
 dup.pair[i,2] <- sample.ID[dup.cordinate[i,2]]
}

rem <- unique(dup.pair[,2])

# save data
save(dup.pair, file="dup.pair.Rdata")
save(rem, file="rem.Rdata")



