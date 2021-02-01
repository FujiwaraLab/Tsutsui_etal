#Making heatmap of the FPKM expression patterns of ECM genes in the isolated HF populations

# Installation of the packages

library(magrittr)
library(DESeq2)
library(EDASeq)
library(biomaRt)
library(gplots)

# Read ECMM gene raw count data
ECMgene <- read.csv("stringtie_ECMgene_count_matrix_forR.csv", header = TRUE, row.names = 1)

#get gene length information
ECMgenelist <- rownames(ECMgene)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ECMgeneensembl <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "mgi_symbol", values = ECMgenelist, mart = ensembl)
genelength <- getGeneLengthAndGCContent(ECMgeneensembl$ensembl_gene_id, "mmusculus", mode = "biomart")
genelength <- as.data.frame(genelength)
nrow(genelength)
#only 274 genes were returned.
#which exist in ECMgenelist but not in ECMgeneensembl?
setdiff(ECMgenelist, ECMgeneensembl$external_gene_name)
#substitution
ECMgenelist <- sub("5430419D17Rik", "Cdcp3", ECMgenelist)
ECMgenelist <- sub("Ctgf", "Ccn2", ECMgenelist)
ECMgenelist <- sub("Cyr61", "Ccn1", ECMgenelist)
ECMgenelist <- sub("Ddx26b", "Ints6l", ECMgenelist)
ECMgenelist <- sub("Nov", "Ccn3", ECMgenelist)
ECMgenelist <- sub("Vwa9", "Ints14", ECMgenelist)
ECMgenelist <- sub("Wisp1", "Ccn4", ECMgenelist)
ECMgenelist <- sub("Wisp2", "Ccn5", ECMgenelist)
ECMgenelist <- sub("Wisp3", "Ccn6", ECMgenelist)
ECMgeneensembl <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "mgi_symbol", values = ECMgenelist, mart = ensembl)
nrow(ECMgeneensembl)
nrow(ECMgene)
#still different number but now 281 < 283?
Wgenes <- duplicated(ECMgeneensembl["external_gene_name"])
Wgenes
which(Wgenes)
#Col4a6 and Cthrc1 are duplicated!
ECMgeneensembl2 <- ECMgeneensembl[-63,]
ECMgeneensembl3 <- ECMgeneensembl2[-86,]
ECMgeneensembl3
genelength <- getGeneLengthAndGCContent(ECMgeneensembl3$ensembl_gene_id, "mmusculus", mode = "biomart")
nrow(genelength)
genelength <- as.data.frame(genelength)

#set experimental condition
condition <- c("Basal", "Basal", "Basal", "LowerIsthmus", "LowerIsthmus", "LowerIsthmus", "UpperBulge", 
               "UpperBulge", "UpperBulge", "MidBulge", "MidBulge", "MidBulge", "HairGerm", "HairGerm", 
               "HairGerm", "HairGerm", "DermalPapilla", "DermalPapilla", "DermalPapilla", "DermalPapilla", 
               "panDermalFibroblast", "panDermalFibroblast", "panDermalFibroblast")
tissuetype <- c("Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", 
                "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", "Epidermis", 
                "Epidermis", "Epidermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis", "Dermis")
sample_info <- data.frame(condition, tissuetype, row.names = names(ECMgene))
sample_info

#generate the DESeqDataSet
ECMDESeq.ds <- DESeqDataSetFromMatrix(countData = ECMgene,
                                      colData = sample_info,
                                      design = ~ condition)
colData(ECMDESeq.ds) %>% head
assay(ECMDESeq.ds, "counts") %>% head
rowData(ECMDESeq.ds) %>% head
counts(ECMDESeq.ds) %>% str
mcols(ECMDESeq.ds)$basepairs <- genelength$length

#calculate FPKM and normalized count
ECMFPKM <- fpkm(ECMDESeq.ds, robust = TRUE)
head(ECMFPKM)
logFPKM <- log(ECMFPKM+0.01, 2)
write.csv(logFPKM, "DESeq2ECMlog2FPKM.csv")

#upload sorted ECM gene data
orderedECMlogFPKM <- read.csv("geneDESeq2logfpkmintegratedECM.csv", header = TRUE, row.names = 1)
head(orderedECMlogFPKM)

#heatmap representation
dev.new()
heatmap.2(as.matrix(orderedECMlogFPKM), trace = "none", Colv = NULL, Rowv = NULL,
          scale = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = F, density.info = "none",
          cexRow = 0.1, lwid = c(1,4), lhei = c(1,5), keysize = 1, key.title = NA,
          col = colorRampPalette(c("white","yellow","orange", "red", "dark red"))(256))
dev.off()

logFPKM <- log(ECMFPKM+0.1, 2)
write.csv(logFPKM, "DESeq2ECMlog0.1FPKM.csv")
#upload sorted ECM gene data
orderedECMlog0.1FPKM <- read.csv("geneDESeq2log0.1fpkmintegratedECM.csv", header = TRUE, row.names = 1)
dev.new()
heatmap.2(as.matrix(orderedECMlog0.1FPKM), trace = "none", Colv = NULL, Rowv = NULL,
          scale = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = F, density.info = "none",
          cexRow = 0.1, lwid = c(1,3), lhei = c(1,3), keysize = 1, key.title = NA,
          col = colorRampPalette(c("white","yellow","orange", "red", "dark red"))(256))
