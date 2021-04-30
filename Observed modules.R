#!/usr/bin/env Rscript
inputargs = commandArgs(trailingOnly=TRUE)
workingData <- inputargs[2]

#Across cell-type WGCNA, PITT tetrad cohort
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
load("PITT tetrad differential WGCNA input (normalized counts, no trans, no quantile norm).rData")

#Check ExpData's and FeatureData(clinicalData)
#PV
sum(rownames(PVMDDexpData) %in% rownames(PVclinicalData))==nrow(PVMDDexpData)
sum(rownames(PVBipolarexpData) %in% rownames(PVclinicalData))==nrow(PVBipolarexpData)
sum(rownames(PVSCHIZexpData) %in% rownames(PVclinicalData))==nrow(PVSCHIZexpData)
sum(rownames(PVControlexpData) %in% rownames(PVczlinicalData))==nrow(PVControlexpData)

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

#Use beta of 6 for all networks, deepsplit=2
prefix <- gsub("expData", "", workingData)

assign(paste0(prefix, "_network"),
       blockwiseModules(get(workingData), power = 6,
                        networkType = "unsigned", 
                        corType="bicor",
                        maxPOutliers = 0.10,
                        TOMType = "unsigned",
                        deepSplit = 2,
                        minCoreKME = 0.3,
                        minKMEtoStay = 0.1,
                        numericLabels = TRUE, 
                        pamRespectsDendro = TRUE,
                        verbose = 3, 
                        saveTOMs = TRUE,
                        saveTOMFileBase = prefix,
                        maxBlockSize=30000))

#Produce module-level outputs
workingNet <- get(paste0(prefix, "_network"))

assign(paste0(prefix, "_moduleLabels"),
       workingNet$colors)

assign(paste0(prefix, "_moduleColors"),
       labels2colors(workingNet$colors))

assign(paste0(prefix, "_MEs"),
       workingNet$MEs)

assign(paste0(prefix, "_geneTree"),
       workingNet$dendrograms[[1]])

save(list=c(paste0(prefix, "_moduleLabels"), paste0(prefix, "_moduleColors"), paste0(prefix, "_MEs"), paste0(prefix, "_geneTree")), file=paste0(prefix,"_network_construction_auto.rData"))
       