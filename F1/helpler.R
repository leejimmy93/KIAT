# Title: helpler file 
# ========================================================
#   
# This is a helpler file that stores all the function used for F1 data anlysis 

library(vcfR)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library("rtracklayer")
library("tidyverse")

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
  
  temp2 <- data.frame(cbind(vcfchrom,vcfpos,vcfref,vcfalt,vcfgts,vcfgq,vcfdp,vcfro,vcfao))
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

#SNP annotation, get gene ID for SNPs with known chromosome and pos information 

SNP.annotation <- function(SNP_data){
  ### get gff file with gene chrom & pos info, gff3 file must be sorted 
  gff.mRNA <- read.table("~/Reference/B.napus/gff.mRNA")
  colnames(gff.mRNA) <- c("CHROM", "start", "end", "name") 
  
  genes <- GRanges(seqnames = Rle(gff.mRNA$CHROM),
                   ranges = IRanges(start = gff.mRNA$start, end = gff.mRNA$end), 
                   names = gff.mRNA$name)
  
  # SNP information 
  SNP <- GRanges(seqnames = Rle(SNP_data$CHROM), 
                 ranges = IRanges(start = SNP_data$POS, SNP_data$POS), 
                 ID = paste(SNP_data$CHROM, SNP_data$POS, sep = "_"))
  
  # overlap SNP position with gene range 
  SNP_gene <- mergeByOverlaps(SNP, genes)
  SNP_gene_df <- as.data.frame(SNP_gene)
  
  # reform output 
  SNP_gene_df2 <- SNP_gene_df[,c("SNP.seqnames", "SNP.start", "genes.names")] 
  colnames(SNP_gene_df2) <- c("CHROM", "POS", "gene_ID")
  
  SNP_gene_final <- 
  SNP_data %>% 
    left_join(SNP_gene_df2, by=c("CHROM", "POS")) 
  
  return(SNP_gene_final)  
}

## binomial test, add binomial test result (p-value, <0.05 means deviated from binomial distribution) as extra columns 
binomial_test <- function(F1){
  F1.414 <- 
    lapply(1:nrow(F1), function(SNP){
      stats <- binom.test(F1$F1_414_RO[SNP], (F1$F1_414_AO[SNP]+F1$F1_414_RO[SNP]), 0.5)
      pvalue <- stats$p.value
      return(pvalue)
    }
    )

  F1.415 <- 
    lapply(1:nrow(F1), function(SNP){
      stats <- binom.test(F1$F1_415_RO[SNP], (F1$F1_415_AO[SNP]+F1$F1_415_RO[SNP]), 0.5)
      pvalue <- stats$p.value
      return(pvalue)
    }
    )

  F1$BT_414 <- unlist(F1.414) 
  F1$BT_415 <- unlist(F1.415)
  
  return(F1)
}

# check parent-of-origin-specific expression pattern 
ASE_test <- function(F1){
  for(i in 1:nrow(F1)){
    if(F1$BT_414[i] >= 0.05){ # not ASE expression 
      F1$ASE_414[i] = "N"
    }
    else if(F1$BT_414[i] < 0.05 & 
              (F1$F1_414_RO[i]/(F1$F1_414_AO[i]+F1$F1_414_RO[i]) < 0.5 &
                 F1$Ae_RO[i]/(F1$Ae_AO[i]+F1$Ae_RO[i]) < 0.5) |
              (F1$BT_414[i] < 0.05 &
                 F1$F1_414_RO[i]/(F1$F1_414_AO[i]+F1$F1_414_RO[i]) > 0.5 &
                 F1$Ae_RO[i]/(F1$Ae_AO[i]+F1$Ae_RO[i] > 0.5))) { # same as Ae     
      F1$ASE_414[i] = "Ae"
    }
    else {
      F1$ASE_414[i] = "Ol"
    } 
  } 
  
  for(i in 1:nrow(F1)){
    if(F1$BT_415[i] >= 0.05){ # not ASE expression 
      F1$ASE_415[i] = "N"
    }
    else if(F1$BT_415[i] < 0.05 & 
              (F1$F1_415_RO[i]/(F1$F1_415_AO[i]+F1$F1_415_RO[i]) < 0.5 &
                 F1$Ae_RO[i]/(F1$Ae_AO[i]+F1$Ae_RO[i]) < 0.5) |
              (F1$BT_415[i] < 0.05 &
                 F1$F1_415_RO[i]/(F1$F1_415_AO[i]+F1$F1_415_RO[i]) > 0.5 &
                 F1$Ae_RO[i]/(F1$Ae_AO[i]+F1$Ae_RO[i] > 0.5))) { # same as Ae     
      F1$ASE_415[i] = "Ae"
    }
    else {
      F1$ASE_415[i] = "Ol"
    } 
  } 
  
  return(F1)
}

# add Da-Ae allele ratio 
Ae_ratio_test <- function(F1){
  for(i in 1:nrow(F1)){
    if(F1$Ae_GT[i] == -1){ # Ae same as ref allele 
      F1$Ae_ratio_414[i] = F1$F1_414_RO[i]/(F1$F1_414_RO[i] + F1$F1_414_AO[i])
      F1$Ae_ratio_415[i] = F1$F1_415_RO[i]/(F1$F1_415_RO[i] + F1$F1_415_AO[i])
    }
    else { # Ae same as alt allele 
      F1$Ae_ratio_414[i] = F1$F1_414_AO[i]/(F1$F1_414_RO[i] + F1$F1_414_AO[i])
      F1$Ae_ratio_415[i] = F1$F1_415_AO[i]/(F1$F1_415_RO[i] + F1$F1_415_AO[i])
    } 
  } 
  
  return(F1)
}













