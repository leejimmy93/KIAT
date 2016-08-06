# ORA with GOseq (Brassica napa version)
# prerequisit
library(ShortRead);library(goseq);library(GO.db);library("annotate")
# for ggplot heatmap
library(WGCNA);library(ggplot2);library(reshape2);library(scales)

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
# library (plyr)
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