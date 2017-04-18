### function for BnRNAseq 

####### 1) GO annotaion 
# ORA with GOseq (Brassica napa version)
# prerequisit
library(ShortRead);library(goseq);library(GO.db);library("annotate")
# for ggplot heatmap
library(WGCNA);library(ggplot2);library(reshape2);library(scales); library (plyr)

Bn_cdna<-readDNAStringSet("/Users/ruijuanli/Desktop/Brassica_project/KIAT_RNA_seq/data/Brassica_napus.annotation_v5.cds_modified.fa")
head(Bn_cdna)
bias<-nchar(Bn_cdna)
#names(bias)<-substr(names(Br_cdna),1,9)
names(bias)<-names(Bn_cdna)
length(bias) # 101040

#  bias.data vector must have the same length as DEgenes vector!

# convert to list (onetime procedure)
# Bngo<-read.table("/Users/ruijuanli/Desktop/Brassica_project/KIAT_RNA_seq/data/Brassica_napus_GO",header=FALSE) # read Br_GO.txt in excel and save as csv
# head(Bngo)
# tail(Bngo)
#
# Bngo.list <- tapply(as.character(Bngo$V2),Bngo$V1,c)
# head(Bngo.list)
# save(Bngo.list,file="/Users/ruijuanli/Desktop/Brassica_project/KIAT_RNA_seq/data/Bngo.list.Rdata") # this does not work ub goseq after updating R
# # Using manually entered categories.
# # Calculating the p-values...
# # 'select()' returned 1:1 mapping between keys and columns
# Bngo.DF<-as.data.frame(Bngo.list)
# Bngo.DF$gene<-rownames(Bngo.DF)
# Bngo.DF[1:10,]
# do.call(rbind.data.frame, Bngo.list)
# Bngo.DF2<-do.call(rbind.data.frame,Bngo.list) # ????
# Bngo.DF3 <- ldply (Bngo.list, data.frame)
# names(Bngo.DF3)<-c("gene","GO") #if Brgo.list does not work in goseq, use DF3.


load("/Users/ruijuanli/Desktop/Brassica_project/KIAT_RNA_seq/data/Bngo.list.Rdata")
GOseq.Bn.ORA<-function(genelist,padjust=0.05,ontology="BP") { # return GO enrichment table, padjus, padjust=0.05
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias)
  #print(TF)
  pwf<-nullp(TF,bias.data=bias)
  #print(pwf$DEgenes)
  GO.pval <- goseq(pwf,gene2cat=Bngo.list,use_genes_without_cat=TRUE) # format became different in new goseq version (021111). Does not work (042716)
  #GO.pval <- goseq(pwf,gene2cat=Brgo.DF3,use_genes_without_cat=TRUE) # format became different in new goseq version (021111)

  #head(GO.pval)
  if(ontology=="BP") {
    GO.pval2<-subset(GO.pval,ontology=="BP")
  } else if(ontology=="CC") {
    GO.pval2<-subset(GO.pval,ontology=="CC")
  } else {
    GO.pval2<-subset(GO.pval,ontology=="MF")
  }

  GO.pval2$over_represented_padjust<-p.adjust(GO.pval2$over_represented_pvalue,method="BH")
  if(GO.pval2$over_represented_padjust[1]>padjust) stop("no enriched GO")
  else {
    enriched.GO<-GO.pval2[GO.pval2$over_represented_padjust<padjust,]
    print("enriched.GO is")
    print(enriched.GO)

    ## write Term and Definition
    for(i in 1:dim(enriched.GO)[1]) {
      enriched.GO$Term[i]<-Term(GOTERM[[enriched.GO[i,"category"]]])
      enriched.GO$Definition[i]<-Definition(GOTERM[[enriched.GO[i,"category"]]])
    }
    return(enriched.GO)
  }
}

### 2) expression profile graph drawing w/o annotation
# expression profile graph, need to have voom transformed data
expression.pattern.Bn.parent <- function(ID){
  rownames(ID) <- ID$V1
  data <- as.data.frame(vstMat.parent[(c(which(rownames(vstMat.parent) %in% rownames(ID)))),])
  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")
  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)

  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)

  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))
  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]


  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue", title="")

  p

  return(p)
}

### 3) expression profile drawing w/ annotation
expression.pattern.Bn.parent.with.annot <- function(ID, annotation){
  data <- as.data.frame(vstMat.parent[(c(which(rownames(vstMat.parent) %in% ID))),])
  ## add gene annotation to ID
  rownames(annotation) <- annotation$V1
  data <- merge(data, annotation, by="row.names")

  rownames(data) <- paste(data$Row.names, data$V2, sep = "-")
  data <- data[,-c(1,29:31)]

  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")

  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)
  data.melt.reshape
  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)
  data.melt.reshape

  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))

  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]
  data.melt.reshape

  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue")

  return(p)
}

# expression profile graph, need to have voom transformed data
expression.pattern.Bn.parent <- function(ID){
  rownames(ID) <- ID$V1
  data <- as.data.frame(vstMat.parent[(c(which(rownames(vstMat.parent) %in% rownames(ID)))),])
  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")
  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)

  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)

  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))
  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]


  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue", title="")

  p

  return(p)
}

### 4) expression profile drawing w/ annotation F1
expression.pattern.Bn.F1 <- function(ID){
  rownames(ID) <- ID$V1
  data <- as.data.frame(vstMat.F1[(c(which(rownames(vstMat.F1) %in% rownames(ID)))),])
  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")
  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)

  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)

  # order data to specific orders: young, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","flowering","early-silique","late-silique"))
  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]

  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue", title="")

  p

  return(p)
}

### 5) expression profile drawing w/ annotation
expression.pattern.Bn.F1.with.annot <- function(ID, annotation){
  data <- as.data.frame(vstMat.F1[(c(which(rownames(vstMat.F1) %in% ID))),])
  ## add gene annotation to ID
  rownames(annotation) <- annotation$V1
  data <- merge(data, annotation, by="row.names")

  rownames(data) <- paste(data$Row.names, data$V2, sep = "-")
  data <- data[,-c(1,29:31)]

  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")

  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)
  data.melt.reshape
  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)
  data.melt.reshape

  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))

  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]
  data.melt.reshape

  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue")

  return(p)
}

# expression profile graph, need to have voom transformed data
expression.pattern.Bn.parent <- function(ID){
  rownames(ID) <- ID$V1
  data <- as.data.frame(vstMat.parent[(c(which(rownames(vstMat.parent) %in% rownames(ID)))),])
  data$geneID <- rownames(data)
  data.melt <- melt(data)
  data.melt$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt$variable)
  data.melt$tissue <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt$variable)
  data.melt$rep <-  gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt$variable)
  data.melt$group <- paste(data.melt$gt, data.melt$tissue, data.melt$geneID, sep = "_")

  data.melt.reshape <- reshape(data.melt[,c("value", "rep", "group")], idvar = "group", direction = "wide", timevar = "rep")
  data.melt.reshape$mean <- apply(data.melt.reshape[,c(2:4)], 1, mean, na.rm=TRUE)
  data.melt.reshape$min <- apply(data.melt.reshape[,c(2:4)], 1, max, na.rm=T)
  data.melt.reshape$max <- apply(data.melt.reshape[,c(2:4)], 1, min, na.rm=T)

  data.melt.reshape$gt <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\1",data.melt.reshape$group)
  data.melt.reshape$tissue <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\3",data.melt.reshape$group)
  data.melt.reshape$geneID <- gsub("([[:print:]]+)(_)([[:print:]]+)(_)([[:print:]]+)","\\5",data.melt.reshape$group)

  # order data to specific orders: young, bolting, flowering, early-silique, and late-silique
  data.melt.reshape$tissue <- factor(data.melt.reshape$tissue, levels = c("Young","bolting","flowering","early-silique","late-silique"))
  data.melt.reshape <- data.melt.reshape[order(data.melt.reshape$tissue),]


  p <- ggplot(data = data.melt.reshape)
  p <- p + geom_line(aes(x = factor(tissue), y = mean,group=gt, color=gt))
  p <- p + facet_grid(~geneID~gt)
  p <- p + geom_errorbar(mapping=aes(x=tissue,ymin=min,ymax=max, width=0.25))
  p <- p + theme(axis.text.x=element_text(angle=90),strip.text.y = element_text(angle=0),legend.position="none")
  p <- p + labs(y = "mean expression value", x="tissue", title="")

  p

  return(p)
}

###### process vcf file 
reformat.vcf.F1 <- function(vcf.F1, vcf.header.F1){
  vcf.header.F1 <- sub("#","",vcf.header.F1) #get rid of the pound sign

  vcf.header.F1 <- unlist(strsplit(vcf.header.F1,split="\t"))
  colnames(vcf.F1) <- vcf.header.F1

  # Before splitting add NAs to blank cells
  # Ae
  vcf.F1$Ae[is.na(vcf.F1$Ae)] <- "NA:NA:NA:NA:NA:NA:NA"
  Ae.tmp.unique <- matrix(
    unlist(strsplit(vcf.F1$Ae,split = ":")),
    nrow=nrow(vcf.F1),
    byrow=TRUE
  )
  colnames(Ae.tmp.unique) <- paste("Ae",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  # Ol
  vcf.F1$Ol[is.na(vcf.F1$Ol)] <- "NA:NA:NA:NA:NA:NA:NA"
  Ol.tmp.unique <- matrix(
    unlist(strsplit(vcf.F1$Ol,split = ":")),
    nrow=nrow(vcf.F1),
    byrow=TRUE
  )
  colnames(Ol.tmp.unique) <- paste("Ol",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")

  # 414
  vcf.F1$`414F1`[is.na(vcf.F1$`414F1`)] <- "NA:NA:NA:NA:NA:NA:NA"
  F1_414.tmp.unique <- matrix(
    unlist(strsplit(vcf.F1$`414F1`,split = ":")),
    nrow=nrow(vcf.F1),
    byrow=TRUE
  )
  colnames(F1_414.tmp.unique) <- paste("414F1",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")

  vcf.F1$`415F1`[is.na(vcf.F1$`415F1`)] <- "NA:NA:NA:NA:NA:NA:NA"

  F1_415.tmp.unique <- matrix(
    unlist(strsplit(vcf.F1$`415F1`,split = ":")),
    nrow=nrow(vcf.F1),
    byrow=TRUE
  )
  colnames(F1_415.tmp.unique) <- paste("415F1",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")

  ######
  vcf.F1 <- cbind(vcf.F1,Ae.tmp.unique,Ol.tmp.unique,F1_414.tmp.unique, F1_415.tmp.unique, stringsAsFactors=FALSE)

  vcf.F1[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual",
            "Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual",
            "414F1_tot.depth","414F1_ref.depth","414F1_ref.qual","414F1_alt.depth","414F1_alt.qual",
            "415F1_tot.depth","415F1_ref.depth","415F1_ref.qual","415F1_alt.depth","415F1_alt.qual")] <-
    apply(vcf.F1[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual",
                    "Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual",
                    "414F1_tot.depth","414F1_ref.depth","414F1_ref.qual","414F1_alt.depth","414F1_alt.qual",
                    "415F1_tot.depth","415F1_ref.depth","415F1_ref.qual","415F1_alt.depth","415F1_alt.qual")],
          2,
          as.numeric
    )
  return(vcf.F1)
}

# function to generate basic stats for SNP analysis (GATK version)
### import data and reformat 

SNP.GATK.reformat <- function(vcf, vcf.header){ # input are vcf file from GATK after filtering & header of the vcf file 
  colnames(vcf) <- vcf.header 
  
  # correct the file 
  problem.row <- which(vcf$FORMAT!="GT:AD:DP:GQ:PL")
  vcf.corrected <- vcf[-c(problem.row),]
  
  # Before splitting add NAs to blank cells 
  vcf.corrected$Ae[vcf.corrected$Ae=="./.:.:.:.:."] <- "NA:NA,NA:NA:NA:NA,NA,NA"
  Ae.GATK <- matrix(
    unlist(strsplit(vcf.corrected$Ae,split = ":")), 
    nrow=nrow(vcf.corrected),  
    byrow=TRUE
  ) 
  colnames(Ae.GATK) <- paste("Ae",c("gt","ref.alt.depth","approx.depth","genotype.qual","Phred.score"),sep="_")
  
  vcf.corrected$Ol[vcf.corrected$Ol=="./.:.:.:.:."] <- "NA:NA,NA:NA:NA:NA,NA,NA" 
  Ol.GATK <- matrix(
    unlist(strsplit(vcf.corrected$Ol,split = ":")),
    nrow=nrow(vcf.corrected),  
    byrow=TRUE
  ) 
  colnames(Ol.GATK) <- paste("Ol",c("gt","ref.alt.depth","approx.depth","genotype.qual","Phred.score"),sep="_")
  vcf.reform <- cbind(vcf.corrected, Ae.GATK, Ol.GATK,stringsAsFactors=FALSE) 
  
  vcf.reform[,c("Ae_approx.depth","Ae_genotype.qual",
                              "Ol_approx.depth","Ol_genotype.qual")] <-
    apply(vcf.reform[,c("Ae_approx.depth","Ae_genotype.qual",
                                      "Ol_approx.depth","Ol_genotype.qual")],
          2,
          as.numeric
    )
  return(vcf.reform)  
}   

#### basic filter (based on QUAL score and snpcluster)
SNP.GATK.basic.filter <- function(vcf){
  # filter based on snpcluster 
  vcf.pass <- vcf[vcf$FILTER!="SnpCluster",]
  snpcluster.pass.ratio <- nrow(vcf.pass)/nrow(vcf)
  
  # filter based on QUAL score 
  vcf.HQ <- vcf.pass[vcf.pass$QUAL>40,]
  QUAL.40.pass.ratio <- nrow(vcf.HQ) / nrow(vcf.pass) 
  Ae.gt.matrix <- table(vcf.HQ$Ae_gt)
  Ol.gt.matrix <- table(vcf.HQ$Ol_gt)   
  
  cat("the percentage of SNPs that are not in snpcluster:", snpcluster.pass.ratio, "\n")
  cat("The percentage of SNPs with QUAL > 40:", QUAL.40.pass.ratio, "\n")
  cat("genotyping call matrix for Ae:", "\n")
  print(Ae.gt.matrix)
  cat("genotyping call matrix for Ol:", "\n")
  print(Ol.gt.matrix)
  
  return(vcf.HQ)
}

###### # function to generate basic stats for SNP analysis (freebayes version)
### import data and reformat 
SNP.freebayes.reformat <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)

  vcf$Ae[is.na(vcf$Ae)] <- "NA:NA:NA:NA:NA:NA:NA"

  Ae.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ae,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  )

  colnames(Ae.tmp.unique) <- paste("Ae",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")

  vcf$Ol[is.na(vcf$Ol)] <- "NA:NA:NA:NA:NA:NA:NA"

    Ol.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ol,split = ":")),
    nrow=nrow(vcf),
    byrow = TRUE
  )

  colnames(Ol.tmp.unique) <- paste("Ol",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")

  vcf.reform <- cbind(vcf,Ae.tmp.unique,Ol.tmp.unique,stringsAsFactors=FALSE)

  vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual","Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")] <- 
    apply(vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual","Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")],
          2,
          as.numeric) 

    return(vcf.reform)  
}

###### # function to generate basic stats for SNP analysis (freebayes version for single sample)
### import data and reformat 
SNP.freebayes.reformat.Ae <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)
  
  vcf$Ae[is.na(vcf$Ae)] <- "NA:NA:NA:NA:NA:NA:NA"
  
  Ae.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ae,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  )
  
  colnames(Ae.tmp.unique) <- paste("Ae",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  
  vcf.reform <- cbind(vcf,Ae.tmp.unique,stringsAsFactors=FALSE)
  
  vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")] <- 
    apply(vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")],
          2,
          as.numeric) 
  
  return(vcf.reform)  
}

SNP.freebayes.reformat.Ol <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)
  
  vcf$Ol[is.na(vcf$Ol)] <- "NA:NA:NA:NA:NA:NA:NA"
  
  Ol.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ol,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  )
  
  colnames(Ol.tmp.unique) <- paste("Ol",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  
  vcf.reform <- cbind(vcf,Ol.tmp.unique,stringsAsFactors=FALSE)
  
  vcf.reform[,c("Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")] <- 
    apply(vcf.reform[,c("Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")],
          2,
          as.numeric) 
  
  return(vcf.reform)  
}

### with GQ 
SNP.freebayes.reformat.Ae.GQ <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)
  
  vcf$unknown[is.na(vcf$unknown)] <- "NA:NA:NA:NA:NA:NA:NA:NA:NA"
  
  Ae.tmp.unique <- matrix(
    unlist(strsplit(vcf$unknown,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  )
  
  
  
  colnames(Ae.tmp.unique) <- paste("Ae",c("gt", "gt.qual","tot.depth", "allele.obs","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  
  vcf.reform <- cbind(vcf,Ae.tmp.unique,stringsAsFactors=FALSE)
  
  vcf.reform[,c("Ae_gt.qual","Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")] <- 
    apply(vcf.reform[,c("Ae_gt.qual","Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")],
          2,
          as.numeric) 
  
  return(vcf.reform)  
}

SNP.freebayes.reformat.Ol <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)
  
  vcf$Ol[is.na(vcf$Ol)] <- "NA:NA:NA:NA:NA:NA:NA"
  
  Ol.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ol,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  )
  
  colnames(Ol.tmp.unique) <- paste("Ol",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  
  vcf.reform <- cbind(vcf,Ol.tmp.unique,stringsAsFactors=FALSE)
  
  vcf.reform[,c("Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")] <- 
    apply(vcf.reform[,c("Ol_tot.depth","Ol_ref.depth","Ol_ref.qual","Ol_alt.depth","Ol_alt.qual")],
          2,
          as.numeric) 
  
  return(vcf.reform)  
} 

##### get fastq sequence for lab test
library(Biostrings) 
SNP.reform <- function(SNP.csv){
  # get SNP data
  SNP.revised <- paste(SNP.csv$CHROM, SNP.csv$POS, sep = "_")
  # get randomly SNPs ... 
  set.seed(100) 
  test <- sample(SNP.revised, size = 100, replace = F) # for 96 plate
  # get 150bp position info flanking the candidate SNPs 
  test.2 <- data.frame(CHROM = gsub("([[:print:]]+)(_)([[:print:]]+)", "\\1", test),
                       POS = gsub("([[:print:]]+)(_)([[:print:]]+)", "\\3", test) 
  ) 
  return(test.2)
  }

get.fasta <- function(test.2, genome.ref.DNAbiostring){
  test.2$start <- as.numeric(as.character(test.2$POS))-150 # 150bp flanking the candidate SNPs
  test.2$end <- as.numeric(as.character(test.2$POS))+150
  # extract 150bp flanking sequence in fasta format 
  seq <- list()
  
  for (i in 1:length(genome.ref.DNAbiostring)){
    test.i <- test.2[test.2$CHROM == names(genome.ref.DNAbiostring)[i],]
    seq[[i]] <- DNAStringSet(napus[[i]], start = test.i$start, end = test.i$end, use.names = T) 
    names(seq[[i]]) <- paste(test.i$CHROM, test.i$POS, sep = "_")
  }
  # merge DNAstringset into a large one 
  seq.final <- do.call("c", seq) 
  return(seq.final) 
} 

#### get GATK vcf content for picked SNPs 
get.vcf <- function(test.2, vcf){
  vcf.new <- vcf[((vcf$CHROM %in% test.2$CHROM) & (vcf$POS %in% test.2$POS)),]
  return(vcf.new) 
}
 






























