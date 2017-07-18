# Title: helpler file 
# ========================================================
#   
#   This is a helpler file that stores all the function used for F1 data anlysis 

##########################################################
# filter based on depth field
# if the depth for any sample is below n, that site will be filtered out.

DP.filter <- function(vcf, n){
  tmp <- vcf.data[,grep("dp", colnames(vcf.data))]
  vcf.filtered.DP <- vcf[which(apply(tmp, 1, min) > n),]
  
  cat("number of SNPs after DP filter is:")
  cat(dim(vcf.filtered.DP)[1])
  cat("\n")
  
  return(vcf.filtered.DP) 
} 

# filter based on genotype quality field
# if the genotype quality for any sample is below n, that site will be filtered out. 

GQ.filter <- function(vcf, n){
  tmp <- vcf.data[,grep("gen.qual", colnames(vcf.data))]
  vcf.filtered.GQ <- vcf[which(apply(tmp, 1, min) > n),]
  
  cat("number of SNPs after GQ filter is:")
  cat(dim(vcf.filtered.GQ)[1])
  cat("\n")
  
  return(vcf.filtered.GQ) 
} 


##set header
setHeader <- function(vcf){
  vcf.header <- system("grep '#C'" + toString(vcf), intern = TRUE)
  vcf.header <- sub("","",vcf.header)
  
  vcf.header <- unlist(strsplit(vcf.header,split="\t"))
  colnames(vcf.data) <- vcf.header
  
  system("grep '##INFO'" + vcf)
  system("grep '##FORMAT'"+ vcf)
  
  return(vcf.data)
}

##split and set headers
split <- function(vcf, x){
  vcf$x[is.na(vcf$x)] <-"NA:NA:NA:NA:NA:NA:NA:NA"
  tmp <- matrix(
    unlist(strsplit(vcf$x,split = ":")),
    nrow=nrow(vcf),
    byrow=TRUE)

  colnames(tmp) <- paste(x, c("gt","gen.qual","dp","ro","qr","ao","qa","gl"),sep="_")
 return (tmp)
}  
  
##convert columns
#convCol<-function(vcf,name){
#  vcf.data <- cbind(vcf, Ae.tmp, Ol.tmp, F414_early_silique.tmp, F415_early_silique.tmp,stringsAsFactors=FALSE)
#  
#  vcf.data[, paste(Ae,name,sep="_") paste(Ol,name,sep="_") paste(F414,name,sep="_") paste(F1415,name,sep="_")] <- 
#    apply(vcf.data[, paste(Ae,name,sep="_") paste(Ol,name,sep="_") paste(F414,name,sep="_") paste(F1415,name,sep="_")],
#          2,
#          as.numeric
#    )
#  return vcf.data
#}

# extract all the SNPs from the filtered dataset (vcf.data.filter.DP) 
#to get loci that are homozygous for parents, see whether they are 
#heterozygous in both F1s. 

GT.filter <- function(vcf){
  vcf.GT <- subset(vcf, (((Ae_gt=="0/0" & Ol_gt=="1/1")) | ((Ae_gt=="1/1" & Ol_gt=="0/0"))))
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
