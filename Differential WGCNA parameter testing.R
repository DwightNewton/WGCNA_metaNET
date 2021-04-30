#!/usr/bin/env Rscript
inputargs = commandArgs(trailingOnly=TRUE)

#Differential WGCNA parameter testing
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
#As with consensus, re-save file to be compatible with R/3.4.3 (openBLAS cluster version for quick iteration of settings)



#resaving files as V2 compatiable (R version <3.5)
# load("PITT tetrad differential WGCNA input.rData")
# #
# save(PVMDDexpData, PVBipolarexpData, PVSCHIZexpData, PVControlexpData,PYR23MDDexpData, PYR23BipolarexpData, PYR23SCHIZexpData, PYR23ControlexpData,PYR56MDDexpData, PYR56BipolarexpData, PYR56SCHIZexpData, PYR56ControlexpData,SSTMDDexpData, SSTBipolarexpData, SSTSCHIZexpData, SSTControlexpData,VIPMDDexpData, VIPBipolarexpData, VIPSCHIZexpData, VIPControlexpData,  PVclinicalData, PYR23clinicalData, PYR56clinicalData, SSTclinicalData, VIPclinicalData, file="PITT tetrad differential WGCNA input (V2 compatible).rData", version=2)


#If R version >3.5:
load("PITT tetrad differential WGCNA input (V2 compatible).rData")
#Check ExpData's and FeatureData(clinicalData)
#PV
sum(rownames(PVMDDexpData) %in% rownames(PVclinicalData))==nrow(PVMDDexpData)
sum(rownames(PVBipolarexpData) %in% rownames(PVclinicalData))==nrow(PVBipolarexpData)
sum(rownames(PVSCHIZexpData) %in% rownames(PVclinicalData))==nrow(PVSCHIZexpData)
sum(rownames(PVControlexpData) %in% rownames(PVclinicalData))==nrow(PVControlexpData)

#PYR23
sum(rownames(PYR23MDDexpData) %in% rownames(PYR23clinicalData))==nrow(PYR23MDDexpData)
sum(rownames(PYR23BipolarexpData) %in% rownames(PYR23clinicalData))==nrow(PYR23BipolarexpData)
sum(rownames(PYR23SCHIZexpData) %in% rownames(PYR23clinicalData))==nrow(PYR23SCHIZexpData)
sum(rownames(PYR23ControlexpData) %in% rownames(PYR23clinicalData))==nrow(PYR23ControlexpData)

#PYR56
sum(rownames(PYR56MDDexpData) %in% rownames(PYR56clinicalData))==nrow(PYR56MDDexpData)
sum(rownames(PYR56BipolarexpData) %in% rownames(PYR56clinicalData))==nrow(PYR56BipolarexpData)
sum(rownames(PYR56SCHIZexpData) %in% rownames(PYR56clinicalData))==nrow(PYR56SCHIZexpData)
sum(rownames(PYR56ControlexpData) %in% rownames(PYR56clinicalData))==nrow(PYR56ControlexpData)

#SST
sum(rownames(SSTMDDexpData) %in% rownames(SSTclinicalData))==nrow(SSTMDDexpData)
sum(rownames(SSTBipolarexpData) %in% rownames(SSTclinicalData))==nrow(SSTBipolarexpData)
sum(rownames(SSTSCHIZexpData) %in% rownames(SSTclinicalData))==nrow(SSTSCHIZexpData)
sum(rownames(SSTControlexpData) %in% rownames(SSTclinicalData))==nrow(SSTControlexpData)

#VIP
sum(rownames(VIPMDDexpData) %in% rownames(VIPclinicalData))==nrow(VIPMDDexpData)
sum(rownames(VIPBipolarexpData) %in% rownames(VIPclinicalData))==nrow(VIPBipolarexpData)
sum(rownames(VIPSCHIZexpData) %in% rownames(VIPclinicalData))==nrow(VIPSCHIZexpData)
sum(rownames(VIPControlexpData) %in% rownames(VIPclinicalData))==nrow(VIPControlexpData)

#Fix PYR23/PYR56
rownames(PYR23MDDexpData) <- gsub("_", "-", rownames(PYR23MDDexpData))
rownames(PYR23BipolarexpData) <- gsub("_", "-", rownames(PYR23BipolarexpData))
rownames(PYR23SCHIZexpData) <- gsub("_", "-", rownames(PYR23SCHIZexpData))
rownames(PYR23ControlexpData) <- gsub("_", "-", rownames(PYR23ControlexpData))

rownames(PYR56MDDexpData) <- gsub("_", "-", rownames(PYR56MDDexpData))
rownames(PYR56BipolarexpData) <- gsub("_", "-", rownames(PYR56BipolarexpData))
rownames(PYR56SCHIZexpData) <- gsub("_", "-", rownames(PYR56SCHIZexpData))
rownames(PYR56ControlexpData) <- gsub("_", "-", rownames(PYR56ControlexpData))

#Use beta of 5 for all networks - try a number of parameters
setwd("parameter_test")
#Testing
filestart <- gsub("expData", "", inputargs)
resassignThresholdList <- seq(0, 0.1, by=0.025)
mergeCutHeightList <- seq(0.3, 0.5, by= 0.05)
deepSplitList <- c(1:4)


ModuleData <- matrix(ncol = 4, nrow=100)
colnames(ModuleData) <- c("resassignThreshold", "mergeCutHeight", "deepSplit", "grey_size")
a <- 1

for(i in resassignThresholdList){
  for(j in mergeCutHeightList){
    for(k in deepSplitList){
      
      Consensus <- blockwiseModules(get(inputargs), power = 6,
                                    networkType = "unsigned", 
                                    corType="bicor",
                                    TOMType = "unsigned",
                                    minModuleSize = 30,
                                    reassignThreshold = i, 
                                    mergeCutHeight = j,
                                    deepSplit = k,
                                    numericLabels = TRUE, 
                                    pamRespectsDendro = TRUE,
                                    verbose = 3, 
                                    maxBlockSize=30000)
      
      
      ConsensusmergedColors <- labels2colors(Consensus$colors)
      sort(table(ConsensusmergedColors))
      pdf(paste0(filestart, "_differential_reassign(", i, ")_mergeCutHeight(", j, ")_deepSplit(", k, ").pdf"))
      # sizeGrWindow(12, 9)
      # Plot the dendrogram and the module colors underneath
      plotDendroAndColors(Consensus$dendrograms[[1]],ConsensusmergedColors[Consensus$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
      
      ModuleData[a,c(1:3)] <- c(i,j,k)
      ModuleData[a,4] <- table(ConsensusmergedColors)["grey"]
      print(a)
      a <- a+1
    }
  }
}

write.csv(ModuleData, file=paste0(filestart, "_differential_parameterTest.csv"))
