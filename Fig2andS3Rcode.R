#DE BM genes between Basal and Pan-dermal fibroblst

#Installation of the packages
library(magrittr)
library(DESeq2)
library(ggplot2)

#read BM gene raw counts
BMgene <- read.csv("stringtie_BMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)
BM <- BMgene[, c(1:3,21:23)]  #1:3 for Basal and 21:23 for pan-DF
head(BM)

# input cell type names
condition <- c("Basal", "Basal", "Basal", 
               "Pan_dermal_fibroblast", "Pan_dermal_fibroblast", "Pan_dermal_fibroblast")

#set cell type condition
BM_info <- data.frame(condition, row.names = names(BM))
BM_info

#generate the DESeqDataSet from BM gene raw data
BMDESeq.ds <- DESeqDataSetFromMatrix(countData = BM, colData = BM_info, design = ~condition)

#remove genes without any counts
BMDESeq.ds <- BMDESeq.ds[ rowSums(counts(BMDESeq.ds)) > 0, ]

#calculate the size factor and add it to the data set
BMDESeq.ds <- estimateSizeFactors(BMDESeq.ds)
sizeFactors(BMDESeq.ds)

#estimate dispersions
BMDESeq.ds <- estimateDispersions(BMDESeq.ds)

#count normalization
BMcounts.normalized <- counts(BMDESeq.ds, normalized = TRUE)
head(BMcounts.normalized)
write.csv(BMcounts.normalized, "normalizedBMBasalvsDF.csv")

#calculate DEGs
BMDESeq.ds <- DESeq(BMDESeq.ds)
resultsNames(BMDESeq.ds)

#extract DEG result
BM.res <- results(BMDESeq.ds)
head(BM.res)
BM.sorted.DEG <- BM.res[order(BM.res$padj),]
write.csv(BM.sorted.DEG, "BMDEGresult.csv")

#read ratio data between Basal and DF
BMratio <- read.csv("BMratio2direction.csv", header = TRUE)
BMratio <- data.frame(BMratio)
head(BMratio)

BMratio$order <- factor(BMratio$Gene, as.character(BMratio$Gene))
Basalerror <- aes(ymax = Basal + Basal_SD, ymin = Basal)
DFerror <- aes(ymax = -DF, ymin = -DF - DF_SD)

#set 2 different plot data
gg.Basal <- ggplot(BMratio, aes(order, Basal)) + geom_bar(aes(y = Basal), stat = "identity", fill = "red") +
  geom_errorbar(Basalerror, width = 0.25) + scale_y_continuous(limits = c(0, 300)) + coord_flip()

gg.DF <- ggplot(BMratio, aes(order, DF)) + geom_bar(aes(y = -DF), stat = "identity", fill = "blue") +
  geom_errorbar(DFerror, width = 0.25) + scale_y_continuous(limits = c(-150, 0)) + coord_flip()

dev.new()
print(gg.Basal)
dev.off()

dev.new()
print(gg.DF)
dev.off()

#Stacked bar graph
BMratioS <- read.csv("BMratiostack.csv", header = TRUE)
BMratioS <- data.frame(BMratioS)
BMratioS$Gene <- as.character(BMratioS$Gene)
BMratioS$Gene <- factor(BMratioS$Gene, levels = unique(BMratioS$Gene))
BMratioS$color <- ifelse(BMratioS$Cell == "Basal", "red", "blue")
head(BMratioS)
g <- ggplot(BMratioS, aes(Gene, ratio)) + 
  geom_bar(mapping = aes(x = Gene, y = ratio), stat = "identity",
           fill = BMratioS$color, position = "stack") + coord_flip()
print(g)

#read inter stitial ECM gene raw counts
ISECMgene <- read.csv("stringtie_ISECMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)
ISECM <- ISECMgene[, c(1:3,21:23)]
head(ISECM)

#set cell type condition
ISECM_info <- data.frame(condition, row.names = names(ISECM))
ISECM_info

#generate the DESeqDataSet from BM gene raw data
ISECMDESeq.ds <- DESeqDataSetFromMatrix(countData = ISECM, colData = ISECM_info, design = ~condition)

#remove genes without any counts
ISECMDESeq.ds <- ISECMDESeq.ds[ rowSums(counts(ISECMDESeq.ds)) > 0, ]

#calculate the size factor and add it to the data set
ISECMDESeq.ds <- estimateSizeFactors(ISECMDESeq.ds)
sizeFactors(ISECMDESeq.ds)

#estimate dispersions
ISECMDESeq.ds <- estimateDispersions(ISECMDESeq.ds)

#count normalization
ISECMcounts.normalized <- counts(ISECMDESeq.ds, normalized = TRUE)
head(ISECMcounts.normalized)
write.csv(ISECMcounts.normalized, "normalizedISECMBasalvsDF.csv")

#calculate DEGs
ISECMDESeq.ds <- DESeq(ISECMDESeq.ds)
resultsNames(ISECMDESeq.ds)

#extract DEG result
ISECM.res <- results(ISECMDESeq.ds)
head(ISECM.res)
ISECM.sorted.DEG <- ISECM.res[order(ISECM.res$padj),]
write.csv(ISECM.sorted.DEG, "ISECMDEGresult.csv")

#read ratio data between Basal and DF
ISECMratio <- read.csv("ISECMratio2direction.csv", header = TRUE)
ISECMratio <- data.frame(ISECMratio)
head(ISECMratio)
ISECMratio$order <- factor(ISECMratio$Gene, as.character(ISECMratio$Gene))
Basalerror <- aes(ymax = Basal + Basal_SD, ymin = Basal)
DFerror <- aes(ymax = - DF, ymin = -DF - DF_SD)

#set 2 different plot data
gg.Basal <- ggplot(ISECMratio, aes(order, Basal)) + geom_bar(aes(y = Basal), stat = "identity", fill = "red") +
  geom_errorbar(Basalerror, width = 0.25) + scale_y_continuous(limits = c(0, 150)) + coord_flip()

gg.DF <- ggplot(ISECMratio, aes(order, DF)) + geom_bar(aes(y = -DF), stat = "identity", fill = "blue") +
  geom_errorbar(DFerror, width = 0.25) + scale_y_continuous(limits = c(-150, 0)) + coord_flip()