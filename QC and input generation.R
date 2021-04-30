#### PITT Tetrad WGCNA QC - consensus and CT-specific
library(WGCNA)
options(stringsAsFactors = FALSE)
options(scipen = 999)
allowWGCNAThreads()
library(dplyr)
library(DESeq2)
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(vsn)



setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/")
load("seFullData_proteincoding.rData")
dds0 <- DESeqDataSet(seFullData_pc, design=~Subject.Group)

genemeta <- as.data.frame(rowData(dds0))

#Filter low-expressing genes
dds0PV <- dds0[,grep("PV", colnames(dds0))]
PV_keep <- rowSums(assay(dds0PV)) > 10 & rowSums(assay(dds0PV) == 0) <= 60
dds0PYR23 <- dds0[,grep("2n3", colnames(dds0))]
PYR23_keep <- rowSums(assay(dds0PYR23)) > 10 & rowSums(assay(dds0PYR23) == 0) <= 60
dds0PYR56 <- dds0[,grep("5n6", colnames(dds0))]
PYR56_keep <- rowSums(assay(dds0PYR56)) > 10 & rowSums(assay(dds0PYR56) == 0) <= 60
dds0SST <- dds0[,grep("SST", colnames(dds0))]
SST_keep <- rowSums(assay(dds0SST)) > 10 & rowSums(assay(dds0SST) == 0) <= 60
dds0VIP <- dds0[,grep("VIP", colnames(dds0))]
VIP_keep <- rowSums(assay(dds0VIP)) > 10 & rowSums(assay(dds0VIP) == 0) <= 60

keep_list <- list(PV_keep, PYR23_keep, PYR56_keep, SST_keep, VIP_keep)
AllCT_keep <- Reduce("|", keep_list)

#Get library size (for examined genes)
PVlib <- data.frame("X"=colnames(dds0[,grep("PV", colnames(dds0))]), "libSize"=colSums(assay(dds0[,grep("PV", colnames(dds0))])))
PYR23lib <- data.frame("X"=colnames(dds0[,grep("L2n3", colnames(dds0))]), "libSize"=colSums(assay(dds0[,grep("L2n3", colnames(dds0))])))
PYR56lib <- data.frame("X"=colnames(dds0[,grep("L5n6", colnames(dds0))]), "libSize"=colSums(assay(dds0[,grep("L5n6", colnames(dds0))])))
SSTlib <- data.frame("X"=colnames(dds0[,grep("SST", colnames(dds0))]), "libSize"=colSums(assay(dds0[,grep("SST", colnames(dds0))])))
VIPlib <- data.frame("X"=colnames(dds0[,grep("VIP", colnames(dds0))]), "libSize"=colSums(assay(dds0[,grep("VIP", colnames(dds0))])))


dds0 <- dds0[AllCT_keep,]
dds0 <- estimateSizeFactors(dds0)

# #FPKM
# rowData(dds_PV_exons_groupCT_age_sex_group)
# RD <- rowData(dds_PV_exons_groupCT_age_sex_group)
# RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
# rowData(dds_PV_exons_groupCT_age_sex_group) <- RD
# PV_FPKM <- as.data.frame(fpkm(dds_PV_exons_groupCT_age_sex_group))
# PV_FPKM$Gene.stable.ID <- rownames(PV_FPKM)
# #


ncounts <- counts(dds0, normalized=TRUE)
ncountsData <- as.data.frame(ncounts)

ncountsData <- ncountsData[,-c(380,381)]


rawcounts <- assay(dds0)
sort(colSums(rawcounts))

clinicalData <- as.data.frame(colData(dds0))

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Transdiagnostic Analyses/Consensus WGCNA/Pre-analysis/")

#Consensus first
#PV
PVexpData <- t(ncountsData[,grep("PVALB", names(ncountsData))])
PVclinicalData <- clinicalData[grep("PVALB", rownames(clinicalData)),]
gsgPV <- goodSamplesGenes(PVexpData, verbose=3)

gsgPV

if (!gsgPV$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgPV$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(PVexpData)[!gsgPV$goodGenes], collapse = ", ")));
  if (sum(!gsgPV$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(PVexpData)[!gsgPV$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  PVexpData = PVexpData[gsgPV$goodSamples, gsgPV$goodGenes]
}

#PYR23
PYR23expData <- t(ncountsData[,grep("L2n3", names(ncountsData))])
PYR23clinicalData <- clinicalData[grep("L2n3", rownames(clinicalData)),]
gsgPYR23 <- goodSamplesGenes(PYR23expData, verbose=3)

gsgPYR23

if (!gsgPYR23$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgPYR23$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(PYR23expData)[!gsgPYR23$goodGenes], collapse = ", ")));
  if (sum(!gsgPYR23$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(PYR23expData)[!gsgPYR23$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  PYR23expData = PYR23expData[gsgPYR23$goodSamples, gsgPYR23$goodGenes]
}

#PYR56
PYR56expData <- t(ncountsData[,grep("L5n6", names(ncountsData))])
PYR56clinicalData <- clinicalData[grep("L5n6", rownames(clinicalData)),]
gsgPYR56 <- goodSamplesGenes(PYR56expData, verbose=3)

gsgPYR56

if (!gsgPYR56$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgPYR56$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(PYR56expData)[!gsgPYR56$goodGenes], collapse = ", ")));
  if (sum(!gsgPYR56$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(PYR56expData)[!gsgPYR56$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  PYR56expData = PYR56expData[gsgPYR56$goodSamples, gsgPYR56$goodGenes]
}


#SST
SSTexpData <- t(ncountsData[,grep("SST", names(ncountsData))])
SSTclinicalData <- clinicalData[grep("SST", rownames(clinicalData)),]
gsgSST <- goodSamplesGenes(SSTexpData, verbose=3)

gsgSST

if (!gsgSST$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgSST$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(SSTexpData)[!gsgSST$goodGenes], collapse = ", ")));
  if (sum(!gsgSST$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(SSTexpData)[!gsgSST$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  SSTexpData = SSTexpData[gsgSST$goodSamples, gsgSST$goodGenes]
}


#VIP
VIPexpData <- t(ncountsData[,grep("VIP", names(ncountsData))])
VIPclinicalData <- clinicalData[grep("VIP", rownames(clinicalData)),]
gsgVIP <- goodSamplesGenes(VIPexpData, verbose=3)

gsgVIP

if (!gsgVIP$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgVIP$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(VIPexpData)[!gsgVIP$goodGenes], collapse = ", ")));
  if (sum(!gsgVIP$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(VIPexpData)[!gsgVIP$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  VIPexpData = VIPexpData[gsgVIP$goodSamples, gsgVIP$goodGenes]
}


#Re-check after gene removal
gsgPV <- goodSamplesGenes(PVexpData, verbose=3)
gsgPV
gsgPYR23 <- goodSamplesGenes(PYR23expData, verbose=3)
gsgPYR23
gsgPYR56 <- goodSamplesGenes(PYR56expData, verbose=3)
gsgPYR56
gsgSST <- goodSamplesGenes(SSTexpData, verbose=3)
gsgSST
gsgVIP <- goodSamplesGenes(VIPexpData, verbose=3)
gsgVIP

#Check mean-SD plots after removal
meanSdPlot(t(PVexpData))
meanSdPlot(t(PYR23expData))
meanSdPlot(t(PYR56expData))
meanSdPlot(t(SSTexpData))
meanSdPlot(t(VIPexpData))


datalist <- c("PVexpData", "PYR23expData", "PYR56expData", "SSTexpData", "VIPexpData")
for(i in datalist){print(dim(get(i)))}

for(i in datalist){
  sampleTree <- hclust(dist(as.matrix(get(i))), method="average")
  sizeGrWindow(12,9)
  pdf(paste0("plots_and_sampleClustering_", i,".pdf"))
  par(cex=0.6)
  par(mar=c(0,4,2,0))
  plot(sampleTree, main=paste0("Sample clustering to detect ", i, " outliers"), sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  dev.off()
}

PYR23clinicalData$X <- gsub("_L", "-L", PYR23clinicalData$X)
PYR56clinicalData$X <- gsub("_L", "-L", PYR56clinicalData$X)

#Add library size information
PVclinicalData <- merge(PVclinicalData, PVlib, by="X")
PYR23clinicalData <- merge(PYR23clinicalData, PYR23lib, by="X")
PYR56clinicalData <- merge(PYR56clinicalData, PYR56lib, by="X")
SSTclinicalData <- merge(SSTclinicalData, SSTlib, by="X")
VIPclinicalData <- merge(VIPclinicalData, VIPlib, by="X")


traitlist <- c("PVclinicalData", "PYR23clinicalData", "PYR56clinicalData", "SSTclinicalData", "VIPclinicalData")  
for (i in 1:length(datalist)) {
  workingClinicalData <- get(traitlist[i])
  workingClinicalData <- workingClinicalData[,c(2,5,9,10, 12:15, 41,42)]
  workingClinicalData$Subject.Group <- gsub("Control", "0",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("MDD", "1",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("Bipolar", "2",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("SCHIZ", "3",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- as.numeric(workingClinicalData$Subject.Group)
  workingClinicalData$Sex <- gsub("M", "0",workingClinicalData$Sex)
  workingClinicalData$Sex <- gsub("F", "1",workingClinicalData$Sex)
  workingClinicalData$Sex <- as.numeric(workingClinicalData$Sex) 
  
  sampleTree <- hclust(dist(get(datalist[i])), method="average")
  traitColors <- numbers2colors(workingClinicalData, signed=TRUE)
  sizeGrWindow(12,9)
  pdf(paste0("sample_dendrogram_and_traitheatmap_", datalist[i],".pdf"))
  plotDendroAndColors(sampleTree, traitColors, groupLabels=names(workingClinicalData), main="Sample dendrogram and trait heatmap")
  dev.off()
  
  #Save dendrogram order and do spearman correlations with traits of interest
  orderDF <- data.frame("X"=rownames(get(datalist[i]))[sampleTree$order], "position"=c(1:length(rownames(get(datalist[i]))[sampleTree$order])))
  assign(traitlist[i], merge(get(traitlist[i]), orderDF, by="X"))
}

#Spearman correlations
library(ggplot2)
library(psych)
library(gridExtra)
for(i in 1:length(traitlist)){
  workingClinicalData <- get(traitlist[i])
  print(traitlist[i])
  workingClinicalData <- workingClinicalData[,c(2,5,9,10, 12:15, 41:43)]
  workingClinicalData$Subject.Group <- gsub("Control", "0",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("MDD", "1",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("Bipolar", "2",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- gsub("SCHIZ", "3",workingClinicalData$Subject.Group)
  workingClinicalData$Subject.Group <- as.numeric(workingClinicalData$Subject.Group)
  workingClinicalData$Sex <- gsub("M", "0",workingClinicalData$Sex)
  workingClinicalData$Sex <- gsub("F", "1",workingClinicalData$Sex)
  workingClinicalData$Sex <- as.numeric(workingClinicalData$Sex) 
  
  print(cor(workingClinicalData, method="spearman"))
  print(corr.test(workingClinicalData, method="spearman"))
  for(j in 1:(ncol(workingClinicalData)-1)){
    plotData <- workingClinicalData[,c("position", names(workingClinicalData)[j])]
    names(plotData) <- c("position", "variable")
    assign(paste0(traitlist[i],"_",names(workingClinicalData)[j]),
            ggplot(plotData, aes(x=position, y=variable))+
            geom_point()+
            ylab(names(workingClinicalData)[j])+
            ggtitle(paste0(traitlist[i],"_",names(workingClinicalData)[j]))+
            geom_smooth(method = "loess"))
    
  }
 
}
PVmerge <- PVclinicalData[,c("HU.", "position")]
PYR23merge <- PYR23clinicalData[,c("HU.", "position")]
PYR56merge <- PYR56clinicalData[,c("HU.", "position")]
SSTmerge <- SSTclinicalData[,c("HU.", "position")]
VIPmerge <- VIPclinicalData[,c("HU.", "position")]

names(PVmerge)[2] <- "PV_position"
names(PYR23merge)[2] <- "PYR23_position"
names(PYR56merge)[2] <- "PYR56_position"
names(SSTmerge)[2] <- "SST_position"
names(VIPmerge)[2] <- "VIP_position"

position_merged <- merge(PVmerge, PYR23merge, by="HU.")
position_merged <- merge(position_merged, PYR56merge, by="HU.")
position_merged <- merge(position_merged, SSTmerge, by="HU.")
position_merged <- merge(position_merged, VIPmerge, by="HU.")

cor(position_merged, method="spearman")

#Layout plots
png("PV_cluster_order_correlations.png", width=1200, height=600)
grid.arrange(PVclinicalData_Tetrad, PVclinicalData_Subject.Group, PVclinicalData_Sex, PVclinicalData_Age, PVclinicalData_PMI, PVclinicalData_RIN, PVclinicalData_sizeFactor, PVclinicalData_libSize, ncol=4)
dev.off()

png("PYR23_cluster_order_correlations.png", width=1200, height=600)
grid.arrange(PYR23clinicalData_Tetrad, PYR23clinicalData_Subject.Group, PYR23clinicalData_Sex, PYR23clinicalData_Age, PYR23clinicalData_PMI, PYR23clinicalData_RIN, PYR23clinicalData_sizeFactor, PYR23clinicalData_libSize, ncol=4)
dev.off()

png("PYR56_cluster_order_correlations.png", width=1200, height=600)
grid.arrange(PYR56clinicalData_Tetrad, PYR56clinicalData_Subject.Group, PYR56clinicalData_Sex, PYR56clinicalData_Age, PYR56clinicalData_PMI, PYR56clinicalData_RIN, PYR56clinicalData_sizeFactor, PYR56clinicalData_libSize, ncol=4)
dev.off()

png("SST_cluster_order_correlations.png", width=1200, height=600)
grid.arrange(SSTclinicalData_Tetrad, SSTclinicalData_Subject.Group, SSTclinicalData_Sex, SSTclinicalData_Age, SSTclinicalData_PMI, SSTclinicalData_RIN, SSTclinicalData_sizeFactor, SSTclinicalData_libSize, ncol=4)
dev.off()

png("VIP_cluster_order_correlations.png", width=1200, height=600)
grid.arrange(VIPclinicalData_Tetrad, VIPclinicalData_Subject.Group, VIPclinicalData_Sex, VIPclinicalData_Age, VIPclinicalData_PMI, VIPclinicalData_RIN, VIPclinicalData_sizeFactor, VIPclinicalData_libSize, ncol=4)
dev.off()

###Removing subject Hu1044 as they were a strong outlier in PYR56-cells, and soft outliers/showing a trend in the same direction in all other cell-types
PVexpData <- PVexpData[!grepl("Hu1044\\.bam", rownames(PVexpData)),]
PYR23expData <- PYR23expData[!grepl("Hu1044\\.bam", rownames(PYR23expData)),]
PYR56expData <- PYR56expData[!grepl("Hu1044\\.bam", rownames(PYR56expData)),]
SSTexpData <- SSTexpData[!grepl("Hu1044\\.bam", rownames(SSTexpData)),]
VIPexpData <- VIPexpData[!grepl("Hu1044\\.bam", rownames(VIPexpData)),]


#################################################### Scale-free topology and mean connectivity plots
powers <- c(1:12)
for (i in datalist){
  sft <- pickSoftThreshold(as.matrix(get(i)), powerVector=powers, verbose=5, networkType="signed",blockSize=2000)
  #Scale-free topology plot
  sizeGrWindow(9,5)
  pdf(paste0("SIGNED_Scalefree_topology_index_", i,".pdf"))
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
       main = paste("Scale independence"));
  
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  dev.off()
  #Mean connectivity plot
  sizeGrWindow(9,5)
  pdf(paste0("Mean_Connectivity_", i,".pdf"))
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}


# beta=5 for all: (scalefree R^2 > 0.9) for all cell-types
# for just ncounts(no log or VST trans): beta = 4 (lowest, or 5 to be consistent)
# Use only genes which are common across all cell-types - 16853 genes
genelist <- intersect(colnames(PVexpData), colnames(PYR23expData))
genelist <- intersect(genelist, colnames(PYR56expData))
genelist <- intersect(genelist, colnames(SSTexpData))
genelist <- intersect(genelist, colnames(VIPexpData))
length(genelist)

PVexpData <- PVexpData[,genelist]
PYR23expData <- PYR23expData[,genelist]
PYR56expData <- PYR56expData[,genelist]
SSTexpData <- SSTexpData[,genelist]
VIPexpData <- VIPexpData[,genelist]

#fix hyphen/underscores in PYR-cell rownames
rownames(PYR23expData) <- gsub("Pyr-", "Pyr_", rownames(PYR23expData))
rownames(PYR23expData)
rownames(PYR56expData) <- gsub("Pyr-", "Pyr_", rownames(PYR56expData))
rownames(PYR56expData)

for(i in datalist){print(dim(get(i)))}
finalsub <- NA
for(i in datalist){
  workingsub <- sapply(rownames(get(i)), gsub, pattern="^.*Hu", replacement="")
  finalsub <- c(finalsub, workingsub)
}

finalsub <- finalsub[-1]

for(i in unique(finalsub)){
  if(length(grep(i, finalsub)) !=5){
    finalsub <- finalsub[!grepl(i, finalsub)]
  }
}
finalsub <- unique(finalsub)
finalsub <- paste0("Hu", finalsub)

for(i in datalist){
  workingrnames <- unique(grep(paste(finalsub, collapse = "|"), rownames(get(i)), value=TRUE))
  assign(i, get(i)[workingrnames,])
}
for(i in datalist){print(dim(get(i)))}
#Recheck
gsgPV <- goodSamplesGenes(PVexpData, verbose=3)
gsgPV
gsgPYR23 <- goodSamplesGenes(PYR23expData, verbose=3)
gsgPYR23
gsgPYR56 <- goodSamplesGenes(PYR56expData, verbose=3)
gsgPYR56
gsgSST <- goodSamplesGenes(SSTexpData, verbose=3)
gsgSST
gsgVIP <- goodSamplesGenes(VIPexpData, verbose=3)
gsgVIP



save(PVexpData, PYR23expData, PYR56expData, SSTexpData, VIPexpData, PVclinicalData, PYR23clinicalData, PYR56clinicalData, SSTclinicalData, VIPclinicalData, file="Consensus WGCNA input files(normalized).rData", version=2)
