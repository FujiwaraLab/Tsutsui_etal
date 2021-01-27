# To calculate correlation coefficients among gene expression patterns of isolated populations using all genes, matrisome genes, non-matrisome genes, BM genes and interstitial ECM genes

# Installation of packages

ibrary(magrittr)
library(DESeq2)
library(vsn)
library(ggplot2)
library(corrplot)

# Load data file
# 1. all genes

allgene <- read.csv("stringtie_gene_count_matrix_forR.csv", header = TRUE, row.names = 1)
head(allgene)

# Set experimental conditions

condition <- c("Basal", "Basal", "Basal", "LowerIsthmus", "LowerIsthmus", "LowerIsthmus", "UpperBulge", 
               "UpperBulge", "UpperBulge", "MidBulge", "MidBulge", "MidBulge", "HairGerm", "HairGerm", 
               "HairGerm", "HairGerm", "DermalPapilla", "DermalPapilla", "DermalPapilla", "DermalPapilla", 
               "panDermalFibroblast", "panDermalFibroblast", "panDermalFibroblast")
tissuetype <- c("Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", 
                "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", 
                "Epidermis", "Epidermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis")

sample_info <- data.frame(condition, tissuetype, row.names = names(allgene))
sample_info

# Generate the DESeqDataSet
allDESeq.ds <- DESeqDataSetFromMatrix(countData = allgene,
                                      colData = sample_info,
                                      design = ~ condition)
colData(allDESeq.ds) %>% head
assay(allDESeq.ds, "counts") %>% head
rowData(allDESeq.ds) %>% head
counts(allDESeq.ds) %>% str

# Remove genes without any counts
allDESeq.ds <- allDESeq.ds[ rowSums(counts(allDESeq.ds)) > 0, ]
colSums(allgene)
colSums(counts(allDESeq.ds))

# calculate the size factor and add it to the data set
allDESeq.ds <- estimateSizeFactors(allDESeq.ds)
sizeFactors(allDESeq.ds)
colData(allDESeq.ds)

# get (logged-) normalized counts
allcounts.sf_normalized <- counts(allDESeq.ds, normalized = TRUE)
alllog.norm.counts <- log2(allcounts.sf_normalized + 1)

# visualization of normalized all gene data
par(mfrow=c(2,1))
boxplot(allcounts.sf_normalized, notch = TRUE, main = "untransformed read counts", ylab = "read counts")
boxplot(alllog.norm.counts, notch = TRUE, 
        main = "log2-transformed read counts",
        ylab = "log2(read counts)")
plot(alllog.norm.counts[,1:2], cex=.1, main = "Normalized log2(read counts)")

msd_plot <- meanSdPlot(alllog.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot$gg +
  ggtitle("sequencing depth normalized log2(read counts)") +
  ylab("standard deviation")

# rlog normalization of all gene and visualization
allDESeq.rlog <- rlog(allDESeq.ds, blind = TRUE)
allrlog.norm.counts <- assay(allDESeq.rlog)
rlogmsd_plot <- meanSdPlot(allrlog.norm.counts,
                           ranks=FALSE,
                           plot = FALSE)
rlogmsd_plot$gg + ggtitle("rlog-transformed read counts") +
  ylab("staadars deviation")

par(mfrow=c(2,1))
boxplot(alllog.norm.counts, notch = TRUE, 
        main = "log2-transformed read counts",
        ylab = "log2(read counts)")
boxplot(allrlog.norm.counts, notch = TRUE, 
        main = "rlog-transformed read counts",
        ylab = "rlog(read counts)")

# calculate the correlation between columns for hclust
distanceS.all_rlog <- as.dist(1-cor(allrlog.norm.counts, method = "spearman"))
plot( hclust(distanceS.all_rlog, method = "complete"), hang = -1,
      labels = colnames(allrlog.norm.counts),
      main = "rlog transformed read counts\ndistance: Spearman correlation")
disteu.all_rlog <- dist(t(allrlog.norm.counts), method = "euclidean")
plot( hclust(disteu.all_rlog, method = "complete"), hang = -1,
      labels = colnames(allrlog.norm.counts),
      main = "rlog transformed read counts\ndistance: euclidean")

# PCA plot using prcomp
pcall <- prcomp(t(allrlog.norm.counts))
sall <- summary(pcall)
sall
pch.group <- c(rep(4, times = 3), rep(0, times = 3), rep(2, times = 3), rep(1, times = 3), rep(5, times = 4),
               rep(16, times = 4), rep(17, times = 3))
par(mfrow=c(1,1))
plot(pcall$x[,1], xlab = paste("PC1 (",round(sall$importance[2]*100, 1)," %)", sep = ""),
     pcall$x[,2], ylab = paste("PC2 (",round(sall$importance[5]*100, 1)," %)", sep = ""),
     col = "black", pch= pch.group, cex=2, las=1, asp=1,
     main = "PCA of seq.depth normalized\n and rlog-transformed all gene read counts")
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
legend("topleft", legend = c("Basal", "Lower isthmus", "Upper bulge", "Mid bulge", "Hair germ",
                             "Dermal papilla", "Pan-dermal fibroblast"),
       col = "black", pch=c(4, 0, 2, 1, 5, 16, 17), cex = 0.75)

# PCA plot using DESeq2
allP <- plotPCA(allDESeq.rlog)
allP <- allP + theme_bw() + ggtitle("rlog transformed counts")
print(allP)

# Correlation coefficient analysis
all.cor <- cor(allrlog.norm.counts, method = "spearman")
corrplot(all.cor, method = "color")
write.csv(all.cor, “corallgenes.csv”)
corrplot.mixed(all.cor, upper = "color", lower = "number", tl.cex = 0.2, number.cex = 0.4)

# 2. Same analyses for ECM genes
ECMgene <- read.csv("stringtie_ECMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)
head(ECMgene)

# generate the DESeqDataSet
ECMDESeq.ds <- DESeqDataSetFromMatrix(countData = ECMgene,
                                      colData = sample_info,
                                      design = ~ condition)
colData(ECMDESeq.ds) %>% head
assay(ECMDESeq.ds, "counts") %>% head
rowData(ECMDESeq.ds) %>% head
counts(ECMDESeq.ds) %>% str

# remove genes without any counts
ECMDESeq.ds <- ECMDESeq.ds[ rowSums(counts(ECMDESeq.ds)) > 0, ]
colSums(ECMgene)
colSums(counts(ECMDESeq.ds))

# calculate the size factor and add it to the data set
ECMDESeq.ds <- estimateSizeFactors(ECMDESeq.ds)
sizeFactors(ECMDESeq.ds)
colData(ECMDESeq.ds)

# count normalization by DE.Seq2
ECMcounts.sf_normalized <- counts(ECMDESeq.ds, normalized = TRUE)
write.csv(ECMcounts.sf_normalized, "normlizedDESeq2_ECMcount_matrix.csv")
ECMlog.norm.counts <- log2(ECMcounts.sf_normalized + 1)

# visualization of normalized data
par(mfrow=c(2,1))
boxplot(ECMcounts.sf_normalized, notch = FALSE, main = "untransformed ECM gene read counts", ylab = "read counts")
boxplot(ECMlog.norm.counts, notch = FALSE, 
        main = "log2-transformed ECM gene read counts",
        ylab = "log2(read counts)")
par(mfrow=c(1,1))
plot(ECMlog.norm.counts[,1:2], cex=.1, main = "Normalized log2(read counts)")
msd_plot <- meanSdPlot(ECMlog.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot$gg +
  ggtitle("sequencing depth normalized log2(read counts)") +
  ylab("standard deviation")

# rlog normalization and visualization of ECM gene data
ECMDESeq.rlog <- rlog(ECMDESeq.ds, blind = TRUE)
ECMrlog.norm.counts <- assay(ECMDESeq.rlog)
write.csv(ECMrlog.norm.counts, "generlogECM.csv")

rlogmsd_plot <- meanSdPlot(ECMrlog.norm.counts,
                           ranks=FALSE,
                           plot = FALSE)
rlogmsd_plot$gg + ggtitle("rlog-transformed read counts") +
  ylab("staadars deviation")
par(mfrow=c(2,1))
boxplot(ECMlog.norm.counts, notch = TRUE, 
        main = "log2-transformed ECM gene read counts",
        ylab = "log2(read counts)")
boxplot(ECMrlog.norm.counts, notch = TRUE, 
        main = "rlog-transformed ECM gene read counts",
        ylab = "rlog(read counts)")
par(mfrow=c(1,1))

# calculate the correlation between columns
distanceS.ECM_rlog <- as.dist(1-cor(ECMrlog.norm.counts, method = "spearman"))
plot( hclust(distanceS.ECM_rlog, method = "complete"), hang = -1,
      labels = colnames(ECMrlog.norm.counts),
      main = "rlog transformed ECM gene read counts\ndistance: Spearman correlation")
disteu.ECM_rlog <- dist(t(ECMrlog.norm.counts), method = "euclidean")
plot( hclust(disteu.ECM_rlog, method = "complete"), hang = -1,
      labels = colnames(ECMrlog.norm.counts),
      main = "rlog transformed ECM gene read counts\ndistance: euclidean")

# PCA plot for ECM genes using prcomp
pcECM <- prcomp(t(ECMrlog.norm.counts))
sECM <- summary(pcECM)
sECM

pch.group <- c(rep(4, times = 3), rep(0, times = 3), rep(2, times = 3), rep(1, times = 3), rep(5, times = 4),
               rep(16, times = 4), rep(17, times = 3))
plot(pcECM$x[,1], xlab = paste("PC1 (",round(sECM$importance[2]*100, 1)," %)", sep = ""),
     pcECM$x[,2], ylab = paste("PC2 (",round(sECM$importance[5]*100, 1)," %)", sep = ""),
     col = "black", pch= pch.group, cex=2, las=1, asp=1,
     main = "PCA of seq.depth normalized\n and rlog-transformed ECM gene read counts")
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
legend("topleft", legend = c("Basal", "Lower isthmus", "Upper bulge", "Mid bulge", "Hair germ",
                             "Dermal papilla", "Pan-dermal fibroblast"),
       col = "black", pch=c(4, 0, 2, 1, 5, 16, 17), cex = 0.75)

# PCA plot using DESeq2
ECMP <- plotPCA(ECMDESeq.rlog)
ECMP <- ECMP + theme_bw() + ggtitle("rlog transformed ECM gene counts")
print(ECMP)

# Correlation coefficient analysis
ECM.cor <- cor(ECMrlog.norm.counts, method = "spearman")
corrplot(ECM.cor, method = "color")
write.csv(ECM.cor, “corECMgenes.csv”)
corrplot.mixed(ECM.cor, upper = "color", lower = "number", tl.cex = 0.2, number.cex = 0.4)

# 3. Analysis of non-ECM genes
# making the list of ECM genes to be removed
removegenename <- rownames(ECMgene)
dim(ECMgene)
dim(allgene)

# making data frame without ECM genes
woECMgene <- allgene[!row.names(allgene)%in%removegenename,]
dim(woECMgene)
head(woECMgene)

sample_info <- data.frame(condition, tissuetype, row.names = names(woECMgene))
sample_info

# generate the DESeqDataSet
woECMDESeq.ds <- DESeqDataSetFromMatrix(countData = woECMgene,
                                        colData = sample_info,
                                        design = ~ condition)
colData(woECMDESeq.ds) %>% head
assay(woECMDESeq.ds, "counts") %>% head
rowData(woECMDESeq.ds) %>% head
counts(woECMDESeq.ds) %>% str

# remove genes without any counts
woECMDESeq.ds <- woECMDESeq.ds[ rowSums(counts(woECMDESeq.ds)) > 0, ]

# calculate the size factor and add it to the data set
woECMDESeq.ds <- estimateSizeFactors(woECMDESeq.ds)
sizeFactors(woECMDESeq.ds)
colData(woECMDESeq.ds)

# rlog normalization of Non-ECM gene expressions
woECMDESeq.rlog <- rlog(woECMDESeq.ds, blind = TRUE)
woECMrlog.norm.counts <- assay(woECMDESeq.rlog)

# calculate the correlation between columns
distanceS.woECM_rlog <- as.dist(1-cor(woECMrlog.norm.counts, method = "spearman"))
plot( hclust(distanceS.woECM_rlog, method = "complete"), hang = -1,
      labels = colnames(woECMrlog.norm.counts),
      main = "rlog transformed all but without ECM gene\ndistance: Spearman correlation")
disteu.woECM_rlog <- dist(t(woECMrlog.norm.counts), method = "euclidean")
plot( hclust(disteu.woECM_rlog, method = "complete"), hang = -1,
      labels = colnames(woECMrlog.norm.counts),
      main = "rlog transformed all but without ECM gene\ndistance: euclidean")

# Correlation coefficient analysis
woECM.cor <- cor(woECMrlog.norm.counts, method = "spearman")
corrplot(woECM.cor, method = "color")
write.csv(woECM.cor, “corNonECMgenes.csv”)
corrplot.mixed(woECM.cor, upper = "color", lower = "number", tl.cex = 0.2, number.cex = 0.4)

# PCA plot using prcomp
pcwoECM <- prcomp(t(woECMrlog.norm.counts))
swoECM <- summary(pcwoECM)
swoECM
pch.group <- c(rep(4, times = 3), rep(0, times = 3), rep(2, times = 3), rep(1, times = 3), rep(5, times = 4),
               rep(16, times = 4), rep(17, times = 3))
dev.off()
plot(pcwoECM$x[,1], xlab = paste("PC1 (",round(swoECM$importance[2]*100, 1)," %)", sep = ""),
     pcwoECM$x[,2], ylab = paste("PC2 (",round(swoECM$importance[5]*100, 1)," %)", sep = ""),
     col = "black", pch= pch.group, cex=2, las=1, asp=1,
     main = "PCA of seq.depth normalized\n and rlog-transformed all but without ECM gene read counts")
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
legend("topleft", legend = c("Basal", "Lower isthmus", "Upper bulge", "Mid bulge", "Hair germ",
                             "Dermal papilla", "Pan-dermal fibroblast"),
       col = "black", pch=c(4, 0, 2, 1, 5, 16, 17), cex = 0.7)

# 4. Load rlog counts of BM ECM genes
BMrlog.norm.counts <- read.csv("generlogBM.csv", header = TRUE, row.names = 1)

# calculate the correlation between columns
distanceS.BM_rlog <- as.dist(1-cor(BMrlog.norm.counts, method = "spearman"))
plot( hclust(distanceS.BM_rlog, method = "complete"), hang = -1,
      labels = colnames(BMrlog.norm.counts),
      main = "rlog transformed BM gene read counts\ndistance: Spearman correlation")
disteu.BM_rlog <- dist(t(BMrlog.norm.counts), method = "euclidean")
plot( hclust(disteu.BM_rlog, method = "complete"), hang = -1,
      labels = colnames(BMrlog.norm.counts),
      main = "rlog transformed BM gene read counts\ndistance: euclidean")

# Correlation coefficient analysis
BM.cor <- cor(BMrlog.norm.counts, method = "spearman")
corrplot(BM.cor, method = "color")
write.csv(BM.cor, "corBMgenes.csv", quote = TRUE)

# 5. Load rlog counts of IS ECM genes
ISECMrlog.norm.counts <- read.csv("generlogISECM.csv", header = TRUE, row.names = 1)

# calculate the correlation between columns
distanceS.ISECM_rlog <- as.dist(1-cor(ISECMrlog.norm.counts, method = "spearman"))
plot( hclust(distanceS.ISECM_rlog, method = "complete"), hang = -1,
      labels = colnames(ISECMrlog.norm.counts),
      main = "rlog transformed interstitial ECM gene read counts\ndistance: Spearman correlation")
disteu.ISECM_rlog <- dist(t(ISECMrlog.norm.counts), method = "euclidean")
plot( hclust(disteu.ISECM_rlog, method = "complete"), hang = -1,
      labels = colnames(ISECMrlog.norm.counts),
      main = "rlog transformed interstitial ECM gene read counts\ndistance: euclidean")

# calculate correlation of BM gene expressions  among regions
ISECM.cor <- cor(ISECMrlog.norm.counts, method = "spearman")
corrplot(ISECM.cor, method = "color")
write.csv(ISECM.cor, "corISECMgenes.csv", quote = TRUE)
