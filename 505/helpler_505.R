# functions for 505 analysis 
library(vcfR) 

# reformat vcf file, return numeric values 
reform.vcf <- function(temp){
  vcfbi <- is.biallelic(temp) # return with logics indicating whehter biallelic or not... 
  vcfref <- subset(getREF(temp), subset = vcfbi) # get ref allele
  vcfalt <- subset(getALT(temp), subset = vcfbi) # get alt allele
  vcfchrom <- subset(getCHROM(temp), subset = vcfbi) # get chromosome info 
  vcfpos <- subset(getPOS(temp), subset = vcfbi) # get pos info 
  vcfgts <- subset(extract.gt(temp, element = "GT", IDtoRowNames = F), subset = vcfbi) 
  
  temp2 <- data.frame(cbind(vcfchrom,vcfpos,vcfref,vcfalt,vcfgts))
  colnames(temp2)[1:4] <- c("CHROM","POS","REF","ALT")
  rnames <- rownames(temp2)
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/0","-1",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("0/1","0",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("1/1","1",x)))
  row.names(temp2) <- rnames 
  
  return(temp2) 
}

# reform vcf file, return letter for genotypes 
reform.vcf.atcg <- function(temp){
  vcfbi <- is.biallelic(temp) # return with logics indicating whehter biallelic or not... 
  vcfref <- subset(getREF(temp), subset = vcfbi) # get ref allele
  vcfalt <- subset(getALT(temp), subset = vcfbi) # get alt allele
  vcfchrom <- subset(getCHROM(temp), subset = vcfbi) # get chromosome info 
  vcfpos <- subset(getPOS(temp), subset = vcfbi) # get pos info 
  vcfgts <- subset(extract.gt(temp, element = "GT", IDtoRowNames = F, return.alleles=T), subset = vcfbi) 
  
  temp2 <- data.frame(cbind(vcfchrom,vcfpos,vcfref,vcfalt,vcfgts))
  colnames(temp2)[1:4] <- c("CHROM","POS","REF","ALT")
  rnames <- rownames(temp2)
  row.names(temp2) <- rnames 
  
  return(temp2) 
}

IUPAC.code <- function(temp2){
  temp2 <- data.frame(sapply(temp2, function(x) sub("A/A","A",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("T/T","T",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("C/C","C",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("G/G","G",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("C/T","Y",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("A/G","R",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("A/T","W",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("G/C","S",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("T/G","K",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("C/A","M",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("T/C","Y",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("G/A","R",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("T/A","W",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("C/G","S",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("G/T","K",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("A/C","M",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("\\.","N",x)))
  
  return(temp2)
}

IUPAC.code.reverse <- function(temp2){
  temp2 <- data.frame(sapply(temp2, function(x) sub("A","AA",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("T","TT",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("C","CC",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("G","GG",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("Y","CT",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("R","AG",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("W","AT",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("S","GC",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("K","TG",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("M","CA",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("Y","TC",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("R","GA",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("W","TA",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("S","CG",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("K","GT",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("M","AC",x)))
  temp2 <- data.frame(sapply(temp2, function(x) sub("N","NN",x)))
  
  return(temp2)
}

# this function is for??? 
vcf_sigSNP_extract <- function(vcf, significant_SNP){ 
  vcf.new <-  
    vcf %>% 
    unite(feature2, CHROM, POS)
  
  vcf.new$feature2 <- gsub("chr", "S", vcf.new$feature2)
  vcf.new$feature2 <- toupper(vcf.new$feature2)
  
  vcf_sigSNP <- vcf.new[vcf.new$feature2 %in% significant_SNP,]
  rownames(vcf_sigSNP) <- vcf_sigSNP$feature2
  
  vcf_sigSNP.t <- as.data.frame(t(vcf_sigSNP))
  rownames(vcf_sigSNP.t) <- gsub("\\.","\\-", rownames(vcf_sigSNP.t))
  rownames(vcf_sigSNP.t) <- gsub("X", "", rownames(vcf_sigSNP.t))  
  
  return(vcf_sigSNP.t)
}  

# the function below return the trait gene correlation 

spearman.cor.gene.trait <- function(trait.gene.sub, phenotypes, genes){ 
  spearman.estimate <- data.frame()
  
  for (i in phenotypes){
    for (j in genes){
      spearman.estimate[i,j] <- cor(trait.gene.sub[,j], trait.gene.sub[,i], method="spearman")
    }
  }
  
  # how to get p-value all together 
  spearman.pvalue <- data.frame()
  tmp <- data.frame()
  
  for (i in phenotypes){
    for (j in genes){
      tmp <- cor.test(trait.gene.sub[,j], trait.gene.sub[,i], method="spearman")
      spearman.pvalue[i,j] <- tmp$p.value
    }
  }
  
  spearman.result <- list()
  spearman.result[[1]] <- spearman.estimate # spearman.result[[1]] is cor estimate
  spearman.result[[2]] <- spearman.pvalue # spearman[[2]] is cor p-vlaue 
  
  return(spearman.result)
} 

### pull out expression value in Da-Ae and Da-Ol-1 for important genes 

expression.pattern.Bn.parent.bar <- function(gene, vstMat.parent){
  # rownames(ID) <- gene$V1
  data <- as.data.frame(vstMat.parent[(c(which(rownames(vstMat.parent) %in% gene$V1))),])
  data$geneID <- rownames(data) 
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")
  
  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")
  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  
  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)
  
  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))
  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]
  
  p <- ggplot(data = data.melt.reshape, aes(x = factor(tissue), y = mean, group = gt, fill=gt))
  p <- p + geom_col(position = "dodge")
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max), position=position_dodge(0.9), width=0.25)
  p <- p + facet_wrap(~geneID, nrow = 2) 
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0))
  p <- p + labs(y = "expression value from RNA-seq", x="tissue", title="") 
  
  p  
  
  return(p) 
}  

