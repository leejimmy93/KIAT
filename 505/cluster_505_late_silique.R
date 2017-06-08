# load library


# load data
setwd("")
load("")

# for cluster 
geno.numeric <- data.matrix(temp2.t)
genDist <- as.matrix(dist(geno.numeric))

#perform the multi-dimensional scaling
geno.mds <- as.data.frame(cmdscale(genDist))

plot(geno.mds)
