#read epidermis expression data
ECMdata <- read.csv("epiECMmRNAvsprotein.csv", header = TRUE)
head(ECMdata)

#Spearmans rank correlation analysis
Lamb2 <- ECMdata[ECMdata$Gene == "Lamb2",]
Lamb2
cor(x = Lamb2$mRNA, y = Lamb2$protein, method = "spearman")

Lamb2 <- ECMdata[ECMdata$Gene == "Lamb2",]
Lamb2
cor(x = Lamb2$mRNA, y = Lamb2$protein, method = "pearson")

Rspo3 <- ECMdata[ECMdata$Gene == "Rspo3",]
Rspo3
cor(x = Rspo3$mRNA, y = Rspo3$protein, method = "pearson")

Col4a6 <- ECMdata[ECMdata$Gene == "Col4a6",]
Col4a6
cor(x = Col4a6$mRNA, y = Col4a6$protein, method = "pearson")

Col4a5 <- ECMdata[ECMdata$Gene == "Col4a5",]
Col4a5
cor(x = Col4a5$mRNA, y = Col4a5$protein, method = "pearson")

Crim1 <- ECMdata[ECMdata$Gene == "Crim1",]
Crim1
cor(x = Crim1$mRNA, y = Crim1$protein, method = "pearson")

Lama4 <- ECMdata[ECMdata$Gene == "Lama4",]
Lama4
cor(x = Lama4$mRNA, y = Lama4$protein, method = "pearson")

Tgfbi <- ECMdata[ECMdata$Gene == "Tgfbi",]
Tgfbi
cor(x = Tgfbi$mRNA, y = Tgfbi$protein, method = "pearson")

Col17a1 <- ECMdata[ECMdata$Gene == "Col17a1",]
Col17a1
cor(x = Col17a1$mRNA, y = Col17a1$protein, method = "pearson")

Matn2 <- ECMdata[ECMdata$Gene == "Matn2",]
Matn2
cor(x = Matn2$mRNA, y = Matn2$protein, method = "pearson")

Vwa1 <- ECMdata[ECMdata$Gene == "Vwa1",]
Vwa1
cor(x = Vwa1$mRNA, y = Vwa1$protein, method = "pearson")

Lamc1 <- ECMdata[ECMdata$Gene == "Lamc1",]
Lamc1
cor(x = Lamc1$mRNA, y = Lamc1$protein, method = "pearson")

Papln <- ECMdata[ECMdata$Gene == "Papln",]
Papln
cor(x = Papln$mRNA, y = Papln$protein, method = "pearson")

Col7a1 <- ECMdata[ECMdata$Gene == "Col7a1",]
Col7a1
cor(x = Col7a1$mRNA, y = Col7a1$protein, method = "pearson")

Col4a4 <- ECMdata[ECMdata$Gene == "Col4a4",]
Col4a4
cor(x = Col4a4$mRNA, y = Col4a4$protein, method = "pearson")

Col4a3 <- ECMdata[ECMdata$Gene == "Col4a3",]
Col4a3
cor(x = Col4a3$mRNA, y = Col4a3$protein, method = "pearson")

Aspn <- ECMdata[ECMdata$Gene == "Aspn",]
Aspn
cor(x = Aspn$mRNA, y = Aspn$protein, method = "pearson")

Postn <- ECMdata[ECMdata$Gene == "Postn",]
Postn
cor(x = Postn$mRNA, y = Postn$protein, method = "pearson")

Egfl6 <- ECMdata[ECMdata$Gene == "Egfl6",]
Egfl6
cor(x = Egfl6$mRNA, y = Egfl6$protein, method = "pearson")

Ctgf <- ECMdata[ECMdata$Gene == "Ctgf",]
Ctgf
cor(x = Ctgf$mRNA, y = Ctgf$protein, method = "pearson")

Ltbp3 <- ECMdata[ECMdata$Gene == "Ltbp3",]
Ltbp3
cor(x = Ltbp3$mRNA, y = Ltbp3$protein, method = "pearson")

Igfbp5 <- ECMdata[ECMdata$Gene == "Igfbp5",]
Igfbp5
cor(x = Igfbp5$mRNA, y = Igfbp5$protein, method = "pearson")

Col4a1 <- ECMdata[ECMdata$Gene == "Col4a1",]
Col4a1
cor(x = Col4a1$mRNA, y = Col4a1$protein, method = "pearson")

Hspg2 <- ECMdata[ECMdata$Gene == "Hspg2",]
Hspg2
cor(x = Hspg2$mRNA, y = Hspg2$protein, method = "pearson")

Tnc <- ECMdata[ECMdata$Gene == "Tnc",]
Tnc
cor(x = Tnc$mRNA, y = Tnc$protein, method = "pearson")

Col6a1 <- ECMdata[ECMdata$Gene == "Col6a1",]
Col6a1
cor(x = Col6a1$mRNA, y = Col6a1$protein, method = "pearson")

Col6a2 <- ECMdata[ECMdata$Gene == "Col6a2",]
Col6a2
cor(x = Col6a2$mRNA, y = Col6a2$protein, method = "pearson")

Col4a2 <- ECMdata[ECMdata$Gene == "Col4a2",]
Col4a2
cor(x = Col4a2$mRNA, y = Col4a2$protein, method = "pearson")

Npnt <- ECMdata[ECMdata$Gene == "Npnt",]
Npnt
cor(x = Npnt$mRNA, y = Npnt$protein, method = "pearson")

Col18a1 <- ECMdata[ECMdata$Gene == "Col18a1",]
Col18a1
cor(x = Col18a1$mRNA, y = Col18a1$protein, method = "pearson")

Ltbp2 <- ECMdata[ECMdata$Gene == "Ltbp2",]
Ltbp2
cor(x = Ltbp2$mRNA, y = Ltbp2$protein, method = "pearson")

Smoc1 <- ECMdata[ECMdata$Gene == "Smoc1",]
Smoc1
cor(x = Smoc1$mRNA, y = Smoc1$protein, method = "pearson")

Thsd4 <- ECMdata[ECMdata$Gene == "Thsd4",]
Thsd4
cor(x = Thsd4$mRNA, y = Thsd4$protein, method = "pearson")

Lama2 <- ECMdata[ECMdata$Gene == "Lama2",]
Lama2
cor(x = Lama2$mRNA, y = Lama2$protein, method = "pearson")

Lamc3 <- ECMdata[ECMdata$Gene == "Lamc3",]
Lamc3
cor(x = Lamc3$mRNA, y = Lamc3$protein, method = "pearson")

Col3a1 <- ECMdata[ECMdata$Gene == "Col3a1",]
Col3a1
cor(x = Col3a1$mRNA, y = Col3a1$protein, method = "pearson")

Nid1 <- ECMdata[ECMdata$Gene == "Nid1",]
Nid1
cor(x = Nid1$mRNA, y = Nid1$protein, method = "pearson")

Spon1 <- ECMdata[ECMdata$Gene == "Spon1",]
Spon1
cor(x = Spon1$mRNA, y = Spon1$protein, method = "pearson")

Lama5 <- ECMdata[ECMdata$Gene == "Lama5",]
Lama5
cor(x = Lama5$mRNA, y = Lama5$protein, method = "pearson")

Fn1 <- ECMdata[ECMdata$Gene == "Fn1",]
Fn1
cor(x = Fn1$mRNA, y = Fn1$protein, method = "pearson")

Col6a3 <- ECMdata[ECMdata$Gene == "Col6a3",]
Col6a3
cor(x = Col6a3$mRNA, y = Col6a3$protein, method = "pearson")

Lamb1 <- ECMdata[ECMdata$Gene == "Lamb1",]
Lamb1
cor(x = Lamb1$mRNA, y = Lamb1$protein, method = "pearson")

#read DP expression data
DPdata <- read.csv("DPECMmRNAvsprotein.csv", header = TRUE)
head(DPdata)

Tnn <- DPdata[DPdata$Gene == "Tnn",]
Tnn
cor(x = Tnn$mRNA, y = Tnn$protein, method = "pearson")

Col13a1 <- DPdata[DPdata$Gene == "Col13a1",]
Col13a1
cor(x = Col13a1$mRNA, y = Col13a1$protein, method = "pearson")

Nid2 <- DPdata[DPdata$Gene == "Nid2",]
Nid2
cor(x = Nid2$mRNA, y = Nid2$protein, method = "pearson")

Col6a6 <- DPdata[DPdata$Gene == "Col6a6",]
Col6a6
cor(x = Col6a6$mRNA, y = Col6a6$protein, method = "pearson")

Lama2 <- DPdata[DPdata$Gene == "Lama2",]
Lama2
cor(x = Lama2$mRNA, y = Lama2$protein, method = "pearson")

Lamc3 <- DPdata[DPdata$Gene == "Lamc3",]
Lamc3
cor(x = Lamc3$mRNA, y = Lamc3$protein, method = "pearson")

Spon1 <- DPdata[DPdata$Gene == "Spon1",]
Spon1
cor(x = Spon1$mRNA, y = Spon1$protein, method = "pearson")

Col6a3 <- DPdata[DPdata$Gene == "Col6a3",]
Col6a3
cor(x = Col6a3$mRNA, y = Col6a3$protein, method = "pearson")

Lamb1 <- DPdata[DPdata$Gene == "Lamb1",]
Lamb1
cor(x = Lamb1$mRNA, y = Lamb1$protein, method = "pearson")

Col4a1 <- DPdata[DPdata$Gene == "Col4a1",]
Col4a1
cor(x = Col4a1$mRNA, y = Col4a1$protein, method = "pearson")

Tnc <- DPdata[DPdata$Gene == "Tnc",]
Tnc
cor(x = Tnc$mRNA, y = Tnc$protein, method = "pearson")

Lamc1 <- DPdata[DPdata$Gene == "Lamc1",]
Lamc1
cor(x = Lamc1$mRNA, y = Lamc1$protein, method = "pearson")

Col6a1 <- DPdata[DPdata$Gene == "Col6a1",]
Col6a1
cor(x = Col6a1$mRNA, y = Col6a1$protein, method = "pearson")

Col4a2 <- DPdata[DPdata$Gene == "Col4a2",]
Col4a2
cor(x = Col4a2$mRNA, y = Col4a2$protein, method = "pearson")

Npnt <- DPdata[DPdata$Gene == "Npnt",]
Npnt
cor(x = Npnt$mRNA, y = Npnt$protein, method = "pearson")

Smoc1 <- DPdata[DPdata$Gene == "Smoc1",]
Smoc1
cor(x = Smoc1$mRNA, y = Smoc1$protein, method = "pearson")

Col7a1 <- DPdata[DPdata$Gene == "Col7a1",]
Col7a1
cor(x = Col7a1$mRNA, y = Col7a1$protein, method = "pearson")

Col17a1 <- DPdata[DPdata$Gene == "Col17a1",]
Col17a1
cor(x = Col17a1$mRNA, y = Col17a1$protein, method = "pearson")

Rspo3 <- DPdata[DPdata$Gene == "Rspo3",]
Rspo3
cor(x = Rspo3$mRNA, y = Rspo3$protein, method = "pearson")

Crim1 <- DPdata[DPdata$Gene == "Crim1",]
Crim1
cor(x = Crim1$mRNA, y = Crim1$protein, method = "pearson")
