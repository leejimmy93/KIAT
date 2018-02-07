library(MBASED)
library(tidyverse)

AnnotateSNPs <- function(SNPdata, gff.mRNA){
  # Combines SNP/SNV loci with gene names
  #
  # Args:
  #   SNP.data: SNP data containing positions to be matched with genomic features
  #   gff: A gff file containing only CHROM, START, END, GeneID
  #
  # Returns:
  #   SNP.data with a new column, GeneID
  
  colnames(gff.mRNA) <- c("CHROM", "start", "end", "name") 
  
  genes <- GRanges(seqnames = Rle(gff.mRNA$CHROM),
                   ranges = IRanges(start = gff.mRNA$start, end = gff.mRNA$end), 
                   names = gff.mRNA$name)
  
  SNPs <- GRanges(seqnames = Rle(SNPdata$CHROM), 
                 ranges = IRanges(start = SNPdata$POS, SNPdata$POS), 
                 CHROM = SNPdata$CHROM,
                 POS = SNPdata$POS)
  
  # Overlap SNP position with gene range 
  overlappedGenes <- mergeByOverlaps(SNPs, genes)
  overlappedGenes <- overlappedGenes[, c(2, 3, 5)]
  colnames(overlappedGenes) <- c("CHROM", "POS", "GeneID")
  
  annotatedSNPdata <- SNPdata %>% 
    left_join(as.data.frame(overlappedGenes), by=c("CHROM", "POS")) 
  
  return(annotatedSNPdata)  
}

SummarizeASEResults_1s <- function(MBASEDOutput) {
  # Output: geneOutputDF is an easier way to look at MAF and p-values at the same time
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE = assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[,1])
  
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <-  assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels = unique(lociOutputGR$aseID)))
  return(list(geneOutput=geneOutputDF, locusOutput=lociOutputList))
}

ExtractASE <- function(MBASEDOutput) {
  # Extract only desired genes
  # Modify ASEindexes to vary the strictness of selection.
  # TODO: Implement Benjamini-Hochberg correction (punish p-values based on rank)
  # Currently using Bonferroni as a crude correction measure (punishes all p-values equally)
  
  results <- SummarizeASEResults_1s(MBASEDOutput)

  ASEindexes <- results$geneOutput$pValueASE * 46577 < 0.05 & 
    results$geneOutput$majorAlleleFrequency > 0.7
  
  significantResults <- list(results$geneOutput[ASEindexes, ], 
                             results$locusOutput[ASEindexes, ])
  return(significantResults)
}

SingleSample <- function(annotatedData, mySNVs, genotype){
  # create RangedSummarizedExperiment object as input for runMBASED
  # then runMBASED
  
  RO <- paste(genotype, "RO", sep = "_")
  AO <- paste(genotype, "AO", sep = "_")
  
  mySample <- SummarizedExperiment(
    assays = list(lociAllele1Counts = matrix(annotatedData[, RO], ncol = 1, dimnames = list(names(mySNVs), 'mySample')),
                lociAllele2Counts = matrix(annotatedData[, AO], ncol = 1,  dimnames = list(names(mySNVs), 'mySample'))
                ),
    rowRanges=mySNVs)
  
  MBASEDOutput <- runMBASED(
    ASESummarizedExperiment = mySample,
    numSim = 10000, 
    isPhased = FALSE)  

  return(MBASEDOutput) 
} 

gff.mRNA <- read.table("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/gff.mRNA")

# Data to use
load("/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/mizukikadowaki/project/output/F1.young.GQ.filtered.Rdata")

annotatedData <- AnnotateSNPs(SNPdata = F1.young.GQ.filtered, gff.mRNA = gff.mRNA)
  
# Remove SNVs with no associated genes 
annotatedData <- filter(annotatedData, !is.na(GeneID)) 
mySNVs <- GRanges(
  seqnames = annotatedData$CHROM,
  ranges = IRanges(start = annotatedData$POS, width = 1),
  aseID = as.vector(annotatedData$GeneID),
  allele1 = annotatedData$REF,
  allele2 = annotatedData$ALT)
  
names(mySNVs) <- annotatedData$GeneID 

system.time(
MBASED.F1.414 <- SingleSample(annotatedData, mySNVs, genotype = "F1_414")
)
save(MBASED.F1.414, file = "/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/F1_SNP/ASE/MBASED.F1.414.Rdata")
