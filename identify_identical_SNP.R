# load data
F2_geno_data <- read.table("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/data/F2_Final_SNP_Calls", header = T)
head(F2_geno_data) # print to make sure data is imported 

# reformat data
rownames(F2_geno_data) <- paste(F2_geno_data$CHROM, F2_geno_data$POS, sep = "_")
F2_geno_data <- F2_geno_data[, -c(1:6)]
colnames(F2_geno_data) <- gsub("([[:print:]]+)(\\.)([[:digit:]])", "\\3", colnames(F2_geno_data))

###### remove redundant SNPs
# some SNPs should have exactly the same genotype across all individuals, they are redundant in terms of genetic map construction, so they should be removed. find those SNPs by doing correlation test. 

F2_geno_data_t <- as.data.frame(t(F2_geno_data))
F2_geno_data_t.numeric <- data.matrix(F2_geno_data_t)

# delete markers with all same genotypes across individuals  
test <- apply(F2_geno_data_t.numeric, 2, function(x) length(unique(x[!is.na(x)])))
filter.polymorphsm <- test != 1

F2_geno_data_t.numeric.2 <- F2_geno_data_t.numeric[,filter.polymorphsm]

# correlation test 
options(warn=-1) # suppress warning message
F2_SNP_correlation <- cor(F2_geno_data_t.numeric.2, use = "pairwise.complete.obs")
save(F2_SNP_correlation, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F2/output/F2_SNP_correlation.Rdata") 
options(warn=0) # unsuppress warning message

# print a message 
print("done") 

# calcualte number of SNPs with correlation of 1
nrow(which(F2_SNP_correlation == 1 & lower.tri(F2_SNP_correlation), arr.ind = T, useNames = T))
