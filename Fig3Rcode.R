#for clustering of DESeq2 results

# Installation of the packages

library(magrittr)
library(DESeq2)
library(ggplot2)
library(NMF)

# read a raw data of epidermis
epiECM <- read.csv("stringtie_ECMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)
epiECM <- epiECM[, 1:16]
head(epiECM)

# input cell type names
condition <- c("Basal", "Basal", "Basal", "LowerIsthmus", "LowerIsthmus", "LowerIsthmus", "UpperBulge", 
               "UpperBulge", "UpperBulge", "MidBulge", "MidBulge", "MidBulge", "HairGerm", "HairGerm", 
               "HairGerm", "HairGerm")

#set cell type condition
epi_info <- data.frame(condition, row.names = names(epiECM))
epi_info

#generate the DESeqDataSet from epidermis raw data
epiDESeq.ds <- DESeqDataSetFromMatrix(countData = epiECM, colData = epi_info, design = ~condition)

#remove genes without any counts
epiDESeq.ds <- epiDESeq.ds[ rowSums(counts(epiDESeq.ds)) > 0, ]

#calculate the size factor and add it to the data set
epiDESeq.ds <- estimateSizeFactors(epiDESeq.ds)
sizeFactors(epiDESeq.ds)

#estimate dispersions
epiDESeq.ds <- estimateDispersions(epiDESeq.ds)

#set levels for x-axis order in plotting
epiDESeq.ds$condition <- factor(epiDESeq.ds$condition, levels = c("Basal", "LowerIsthmus", "UpperBulge", "MidBulge", "HairGerm"))

#calculate DEGs using LRT method
epiDESeq.ds <- nbinomLRT(epiDESeq.ds, full = ~condition, reduced = ~ 1)
resultsNames(epiDESeq.ds)

#extract LRT result
res <- results(epiDESeq.ds)
head(res)
epi.sorted.DEG <- res[order(res$padj),]
write.csv(epi.sorted.DEG, "epiLRTDEGresult.csv")

#rlog normalixzation
epiDESeq.rlog <- rlog(epiDESeq.ds, blind = TRUE)
epirlog.norm.counts <- assay(epiDESeq.rlog)
write.csv(epirlog.norm.counts, "generlogepiECM.csv")

#heatmap visualization
epiECMrlog <- read.csv("generlogepiECM.csv", header = TRUE, row.names = 1)
head(epiECMrlog)
epiDEGgenes <- rownames(subset(epi.sorted.DEG, padj < 0.001))
hm.epiDEGgenes <- epiECMrlog[epiDEGgenes,]
dev.new()
aheatmap(hm.epiDEGgenes, Rowv = TRUE, Colv = NA,
         distfun = "spearman", hclustfun = "complete", scale = "row")
dev.off()

#Z-scored rlog data
epiECMrlogZ <- scale(epiECMrlog)
head(epiECMrlogZ)
write.csv(epiECMrlogZ, "epiECMrlogZ.csv")

#load the clustering result containing 10 clustered genes
epiDEGgroup <- read.csv("clusteredgeneinepi.csv", header = TRUE)
head(epiDEGgroup)

#make group 1 data and binarize it by k-means clustering
Zgroup1 <- epiECMrlogZ[epiDEGgroup$group1[1:10],]
Zgroup1
kclu.res.group1Z <- kmeans(t(Zgroup1), centers = 2, iter.max = 100)
kclu.res.group1Z$cluster #Basal, Lower isthmus, HG Vs Upper bulge, Mid bulge

#group 2
Zgroup2 <- epiECMrlogZ[epiDEGgroup$group2[1:10],]
Zgroup2
kclu.res.group2Z <- kmeans(t(Zgroup2), centers = 2, iter.max = 100)
kclu.res.group2Z$cluster #Basal, HG Vs Lower isthmus, Upper bulge, Mid bulge

#group 3
Zgroup3 <- epiECMrlogZ[epiDEGgroup$group3[1:7],]
Zgroup3
kclu.res.group3Z <- kmeans(t(Zgroup3), centers = 2, iter.max = 100)
kclu.res.group3Z$cluster #Upper bulge(Rep2, 3), HG VS Basal, Lower isthmus, Upoer bulge (Rep1), Mid bulge

#group 4
Zgroup4 <- epiECMrlogZ[epiDEGgroup$group4[1:5],]
Zgroup4
kclu.res.group4Z <- kmeans(t(Zgroup4), centers = 2, iter.max = 100)
kclu.res.group4Z$cluster # Lower isthmus, Upper bulge,HG Vs Basal, Mid bulge

#group 5
Zgroup5 <- epiECMrlogZ[epiDEGgroup$group5[1:25],]
Zgroup5
kclu.res.group5Z <- kmeans(t(Zgroup5), centers = 2, iter.max = 100)
kclu.res.group5Z$cluster #Upper bulge, Mid bulge, HG Vs Basal, Lower isthmus

#group 6
Zgroup6 <- epiECMrlogZ[epiDEGgroup$group6[1:32],]
Zgroup6
kclu.res.group6Z <- kmeans(t(Zgroup6), centers = 2, iter.max = 100)
kclu.res.group6Z$cluster #Basal, Lower isthmus, Upper bulge, Mid bulge Vs HG

#group 7
Zgroup7 <- epiECMrlogZ[epiDEGgroup$group7[1:5],]
Zgroup7
kclu.res.group7Z <- kmeans(t(Zgroup7), centers = 2, iter.max = 100)
kclu.res.group7Z$cluster #Basal (Rep1, 3), Lower isthmus, Upper bulge (Rep1, 2) Vs Basal (Rep2), Upper bulge (Rep3), Mid bulge, HG

#group 8
Zgroup8 <- epiECMrlogZ[epiDEGgroup$group8[1:7],]
Zgroup8
kclu.res.group8Z <- kmeans(t(Zgroup8), centers = 2, iter.max = 100)
kclu.res.group8Z$cluster #Upper bulge, Mid bulge Vs Basal, Lower isthmus, HG

#group 9
Zgroup9 <- epiECMrlogZ[epiDEGgroup$group9[1:18],]
Zgroup9
kclu.res.group9Z <- kmeans(t(Zgroup9), centers = 2, iter.max = 100)
kclu.res.group9Z$cluster #Basal, Lower isthmus, Mid bulge (Rep1, 2) Vs Upper bulge, Mid bulge (Rep3), HG

#group 10
Zgroup10 <- epiECMrlogZ[epiDEGgroup$group10[1:9],]
Zgroup10
kclu.res.group10Z <- kmeans(t(Zgroup10), centers = 2, iter.max = 100)
kclu.res.group10Z$cluster #Mid bulge Vs Basal, Lower isthmus, Upper bulge, HG
