# for dermis DEG ECM identification

# Installation of the packages

library(magrittr)
library(DESeq2)
library(ggplot2)
library(NMF)

# read a count data of ECM genes
ECM <- read.csv("stringtie_ECMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)
derECM <- ECM[, 17:23]
head(derECM)

# input cell type names
condition <- c("Dermal_papilla", "Dermal_papilla", "Dermal_papilla", "Dermal_papilla", 
               "Pan_dermal_fibroblast", "Pan_dermal_fibroblast", "Pan_dermal_fibroblast")

#set cell type condition
der_info <- data.frame(condition, row.names = names(derECM))
der_info

#generate the DESeqDataSet from epidermis raw data
derDESeq.ds <- DESeqDataSetFromMatrix(countData = derECM, colData = der_info, design = ~condition)

#remove genes without any counts
derDESeq.ds <- derDESeq.ds[ rowSums(counts(derDESeq.ds)) > 0, ]

#calculate the size factor and add it to the data set
derDESeq.ds <- estimateSizeFactors(derDESeq.ds)
sizeFactors(derDESeq.ds)

#estimate dispersions
derDESeq.ds <- estimateDispersions(derDESeq.ds)

#set levels for x-axis order in plotting
derDESeq.ds$condition <- factor(derDESeq.ds$condition, levels = c("Dermal_papilla", "Pan_dermal_fibroblast"))

#calculate DEGs
derDESeq.ds <- DESeq(derDESeq.ds)
resultsNames(derDESeq.ds)

#extract DEG result
der.res <- results(derDESeq.ds)
head(der.res)
der.sorted.DEG <- der.res[order(der.res$padj),]
write.csv(der.sorted.DEG, "dermisDEGresult.csv")

#rlog normalixzation
derDESeq.rlog <- rlog(derDESeq.ds, blind = TRUE)
derrlog.norm.counts <- assay(derDESeq.rlog)
write.csv(derrlog.norm.counts, "generlogderECM.csv")
#Z-scored rlog data
derECMrlogZ <- scale(derrlog.norm.counts)
write.csv(derECMrlogZ, "derECMrlogZ.csv")

#heatmap visualization
derECMrlog <- read.csv("generlogderECM.csv", header = TRUE, row.names = 1)
head(derECMrlog)
derDEGgenes <- rownames(subset(der.sorted.DEG, padj < 0.001))
hm.derDEGgenes <- derECMrlog[derDEGgenes,]
dev.new()
aheatmap(hm.derDEGgenes, Rowv = TRUE, Colv = NA,
         distfun = "spearman", hclustfun = "complete", scale = "row")
