# Title: helpler file 
# ========================================================
#   
# This is a helpler file that stores all the function used for F1 data anlysis 

library(vcfR)

##########################################################
## formatting vcf file 
# reformat vcf file, return numeric values 
reform.vcf.F1 <- function(temp){
  vcfbi <- is.biallelic(temp) # return with logics indicating whehter biallelic or not... 
  vcfref <- subset(getREF(temp), subset = vcfbi) # get ref allele
  vcfalt <- subset(getALT(temp), subset = vcfbi) # get alt allele
  vcfchrom <- subset(getCHROM(temp), subset = vcfbi) # get chromosome info 
  vcfpos <- subset(getPOS(temp), subset = vcfbi) # get pos info 
  vcfgts <- subset(extract.gt(temp, element = "GT", IDtoRowNames = F), subset = vcfbi)
  vcfgq <- subset(extract.gt(temp, element="GQ", IDtoRowNames=F), subset=vcfbi)
  vcfdp <- subset(extract.gt(temp, element="DP", IDtoRowNames=F), subset=vcfbi) 
  vcfro <- subset(extract.gt(temp, element="RO", IDtoRowNames=F), subset=vcfbi)
  vcfao <- subset(extract.gt(temp, element="AO", IDtoRowNames=F), subset=vcfbi)
  
  temp2 <- data.frame(cbind(vcfchrom,vcfpos,vcfref,vcfalt,vcfgts,vcfgq,vcfdp))
  colnames(temp2)[1:4] <- c("CHROM","POS","REF","ALT")
  colnames(temp2)[5:8] <- paste(c("Ae", "Ol", "F1_414", "F1_415"), c("GT"), sep="_")
  colnames(temp2)[9:12] <- paste(c("Ae", "Ol", "F1_414", "F1_415"), c("GQ"), sep="_")
  colnames(temp2)[13:16] <- paste(c("Ae", "Ol", "F1_414", "F1_415"), c("DP"), sep="_")
  colnames(temp2)[17:20] <- paste(c("Ae", "Ol", "F1_414", "F1_415"), c("RO"), sep="_")
  colnames(temp2)[21:24] <- paste(c("Ae", "Ol", "F1_414", "F1_415"), c("AO"), sep="_")
  
  rnames <- rownames(temp2)
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/0","-1",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/1","0",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("1/1","1",x)))
  row.names(temp2) <- rnames 
  
  return(temp2) 
} 

# filter based on genotype quality field
# if the genotype quality for any sample is below n, that site will be filtered out. 

GQ.filter <- function(vcf, n){
  tmp <- vcf[,grep("GQ", colnames(vcf))]
  vcf.filtered.GQ <- vcf[which(apply(tmp, 1, min) > n),]
  
  cat("number of SNPs after GQ filter is:")
  cat(dim(vcf.filtered.GQ)[1])
  cat("\n")
  
  return(vcf.filtered.GQ) 
} 

# filter based on depth field
# if the depth for any sample is below n, that site will be filtered out.
DP.filter <- function(vcf, n){
  tmp <- vcf[,grep("DP", colnames(vcf))]
  vcf.filtered.DP <- vcf[which(apply(tmp, 1, min) > n),]
  
  cat("number of SNPs after DP filter is:")
  cat(dim(vcf.filtered.DP)[1])
  cat("\n")
  
  return(vcf.filtered.DP) 
} 

# extract all the SNPs from the filtered dataset (vcf.data.filter.DP) 
#to get loci that are homozygous for parents, see whether they are 
#heterozygous in both F1s. 

GT.filter <- function(vcf){
  
  return(vcf.GT)
} 

#SNP annotation 
SNPannotation <- function(gff, vcf){
  library(IRanges)
  library(GenomicRanges)
  library(GenomicFeatures)
  library("rtracklayer") 
  
  gff.mRNA <- read.table(gff)
  colnames(gff.mRNA) <- c("CHROM", "start", "end", "name") 
  
  genes <- GRanges(seqnames = Rle(gff.mRNA$CHROM),ranges = IRanges(start = gff.mRNA$start, end = gff.mRNA$end), names = gff.mRNA$name)
  SNP <- GRanges(seqnames = Rle(vcf$`#CHROM`), ranges = IRanges(start = vcf$POS, end = vcf$POS), ID = paste(vcf$`#CHROM`, vcf$POS, sep = "_"))
  SNP_gene <- mergeByOverlaps(SNP, genes) 
  SNP_gene_df <- as.data.frame(SNP_gene)
  SNP_gene_final <- SNP_gene_df[,c("SNP.ID", "genes.seqnames", "SNP.start", "names")]
  
  return(vcf.filtered.GT)
}
