library(cluster)
library(scatterplot3d)
library(WGCNA)
library(dplyr)
library(igraph)
library(DESeq2)
library(Rsamtools);
library("GenomicFeatures");
library("GenomicAlignments");
library("BiocParallel");
library("pheatmap");
library("RColorBrewer");
library(psych)
options(stringsAsFactors = FALSE)
substrLeft <- function (x, n=5){
  substr(x, nchar(x)-n, nchar(x))}

load("PermutationMasterList.rData")
load("PITT tetrad differential WGCNA input for permutations V2.rData")

cor <- WGCNA::cor

###Outcome measures: #modules (x8), total, negative, and positive connections between meta-networks (CON and UCMS - 12 in total)
MinNum <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
MaxNum <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
group1 <- commandArgs(trailingOnly = TRUE)[3]
group2 <- commandArgs(trailingOnly = TRUE)[4]

ResultsFile <- as.data.frame(matrix(ncol=61, nrow=360))
names(ResultsFile) <- c("PVcon_Mods","PVstr_Mods","PYR23con_Mods","PYR23str_Mods","PYR56con_Mods","PYR56str_Mods","SSTcon_Mods","SSTstr_Mods","VIPcon_Mods","VIPstr_Mods","PVcon_GenesClustered","PVstr_GenesClustered","PYR23con_GenesClustered","PYR23str_GenesClustered","PYR56con_GenesClustered","PYR56str_GenesClustered","SSTcon_GenesClustered","SSTstr_GenesClustered","VIPcon_GenesClustered","VIPstr_GenesClustered","PV_PYR23_Group1_sigconn","PV_PYR56_Group1_sigconn","PV_SST_Group1_sigconn","PV_VIP_Group1_sigconn","PYR23_PYR56_Group1_sigconn","PYR23_SST_Group1_sigconn","PYR23_VIP_Group1_sigconn","PYR56_SST_Group1_sigconn","PYR56_VIP_Group1_sigconn","SST_VIP_Group1_sigconn","PV_PYR23_Group2_sigconn","PV_PYR56_Group2_sigconn","PV_SST_Group2_sigconn","PV_VIP_Group2_sigconn","PYR23_PYR56_Group2_sigconn","PYR23_SST_Group2_sigconn","PYR23_VIP_Group2_sigconn","PYR56_SST_Group2_sigconn","PYR56_VIP_Group2_sigconn","SST_VIP_Group2_sigconn","PV_PYR23_Group1_norm_factor","PV_PYR56_Group1_norm_factor","PV_SST_Group1_norm_factor","PV_VIP_Group1_norm_factor","PYR23_PYR56_Group1_norm_factor","PYR23_SST_Group1_norm_factor","PYR23_VIP_Group1_norm_factor","PYR56_SST_Group1_norm_factor","PYR56_VIP_Group1_norm_factor","SST_VIP_Group1_norm_factor","PV_PYR23_Group2_norm_factor","PV_PYR56_Group2_norm_factor","PV_SST_Group2_norm_factor","PV_VIP_Group2_norm_factor","PYR23_PYR56_Group2_norm_factor","PYR23_SST_Group2_norm_factor","PYR23_VIP_Group2_norm_factor","PYR56_SST_Group2_norm_factor","PYR56_VIP_Group2_norm_factor","SST_VIP_Group2_norm_factor", "density")
setwd(group2)
#Will change to "group 1" and "group 2"
datalist <- ls(pattern="expData")

#Get expression data lists
group1Data <- grep(group1, datalist, value=TRUE)
group2Data <- grep(group2, datalist, value=TRUE)

#Merge the data
PVmerged <- rbind(get(grep("PV", group1Data, value=TRUE)), get(grep("PV", group2Data, value=TRUE)))
PYR23merged <- rbind(get(grep("PYR23", group1Data, value=TRUE)), get(grep("PYR23", group2Data, value=TRUE)))
PYR56merged <- rbind(get(grep("PYR56", group1Data, value=TRUE)), get(grep("PYR56", group2Data, value=TRUE)))
SSTmerged <- rbind(get(grep("SST", group1Data, value=TRUE)), get(grep("SST", group2Data, value=TRUE)))
VIPmerged <- rbind(get(grep("VIP", group1Data, value=TRUE)), get(grep("VIP", group2Data, value=TRUE)))


rowsequence <- c(gsub(".*-+","",rownames(get(group1Data[1]))), gsub(".*-+","",rownames(get(group2Data[1]))))

PVmerged <- PVmerged[order(match(gsub(".*-+","",rownames(PVmerged)), rowsequence)),]
PYR23merged <- PYR23merged[order(match(gsub(".*-+","",rownames(PYR23merged)), rowsequence)),]
PYR56merged <- PYR56merged[order(match(gsub(".*-+","",rownames(PYR56merged)), rowsequence)),]
SSTmerged <- SSTmerged[order(match(gsub(".*-+","",rownames(SSTmerged)), rowsequence)),]
VIPmerged <- VIPmerged[order(match(gsub(".*-+","",rownames(VIPmerged)), rowsequence)),]

gsub(".*-+","",rownames(PVmerged)) == gsub(".*-+","",rownames(PYR23merged))
gsub(".*-+","",rownames(PVmerged)) == gsub(".*-+","",rownames(PYR56merged))
gsub(".*-+","",rownames(PYR56merged)) == gsub(".*-+","",rownames(SSTmerged))
gsub(".*-+","",rownames(SSTmerged)) == gsub(".*-+","",rownames(VIPmerged))

# densitylist <- seq(0.05, 0.4, by=0.01)
start <- seq(0.1, 0.02, by=-0.01)
densitylist <- c(start, start/10, start/100, start/1000)


#Sample sizes are the same across CT's
group1n <- nrow(get(group1Data[1]))
group2n <- nrow(get(group2Data[1]))
NumList <- c(MinNum:MaxNum)

resultsrow <- 1
###change range based on processing chunk
system.time (for(i in 1:length(NumList)){
  
  permutedIDs <- PermutationList[[NumList[i]]]
  
  PVgroup1 <- PVmerged[permutedIDs[1:group1n],]
  PVgroup2 <- PVmerged[permutedIDs[(group1n+1):(group1n+group2n)],]
  PYR23group1 <- PYR23merged[permutedIDs[1:group1n],]
  PYR23group2 <- PYR23merged[permutedIDs[(group1n+1):(group1n+group2n)],]
  PYR56group1 <- PYR56merged[permutedIDs[1:group1n],]
  PYR56group2 <- PYR56merged[permutedIDs[(group1n+1):(group1n+group2n)],]
  SSTgroup1 <- SSTmerged[permutedIDs[1:group1n],]
  SSTgroup2 <- SSTmerged[permutedIDs[(group1n+1):(group1n+group2n)],]
  VIPgroup1 <- VIPmerged[permutedIDs[1:group1n],]
  VIPgroup2 <- VIPmerged[permutedIDs[(group1n+1):(group1n+group2n)],]
  
  PVgroup1net <- blockwiseModules(PVgroup1, 
                                  power = 6,
                                  networkType = "unsigned",
                                  corType="bicor",
                                  maxPOutliers = 0.1,
                                  TOMType = "unsigned", 
                                  minModuleSize = 30,
                                  minCoreKME = 0.3,
                                  minKMEtoStay = 0.1,
                                  deepSplit = 2,
                                  numericLabels = TRUE, 
                                  pamRespectsDendro = TRUE,
                                  saveTOMs = FALSE,
                                  verbose = 1, 
                                  maxBlockSize=3000)
  
  PYR23group1net <- blockwiseModules(PYR23group1, 
                                     power = 6,
                                     networkType = "unsigned",
                                     corType="bicor",
                                     maxPOutliers = 0.1,
                                     TOMType = "unsigned", 
                                     minModuleSize = 30,
                                     minCoreKME = 0.3,
                                     minKMEtoStay = 0.1,
                                     deepSplit = 2,
                                     numericLabels = TRUE, 
                                     pamRespectsDendro = TRUE,
                                     saveTOMs = FALSE,
                                     verbose = 1, 
                                     maxBlockSize=3000)
  
  PYR56group1net <- blockwiseModules(PYR56group1, 
                                     power = 6,
                                     networkType = "unsigned",
                                     corType="bicor",
                                     maxPOutliers = 0.1,
                                     TOMType = "unsigned", 
                                     minModuleSize = 30,
                                     minCoreKME = 0.3,
                                     minKMEtoStay = 0.1,
                                     deepSplit = 2,
                                     numericLabels = TRUE, 
                                     pamRespectsDendro = TRUE,
                                     saveTOMs = FALSE,
                                     verbose = 1, 
                                     maxBlockSize=3000)
  
  SSTgroup1net <- blockwiseModules(SSTgroup1, 
                                   power = 6,
                                   networkType = "unsigned",
                                   corType="bicor",
                                   maxPOutliers = 0.1,
                                   TOMType = "unsigned", 
                                   minModuleSize = 30,
                                   minCoreKME = 0.3,
                                   minKMEtoStay = 0.1,
                                   deepSplit = 2,
                                   numericLabels = TRUE, 
                                   pamRespectsDendro = TRUE,
                                   saveTOMs = FALSE,
                                   verbose = 1, 
                                   maxBlockSize=3000)
  
  VIPgroup1net <- blockwiseModules(VIPgroup1, 
                                   power = 6,
                                   networkType = "unsigned",
                                   corType="bicor",
                                   maxPOutliers = 0.1,
                                   TOMType = "unsigned", 
                                   minModuleSize = 30,
                                   minCoreKME = 0.3,
                                   minKMEtoStay = 0.1,
                                   deepSplit = 2,
                                   numericLabels = TRUE, 
                                   pamRespectsDendro = TRUE,
                                   saveTOMs = FALSE,
                                   verbose = 1, 
                                   maxBlockSize=3000)
  
  
  PVgroup2net <- blockwiseModules(PVgroup2, 
                                  power = 6,
                                  networkType = "unsigned",
                                  corType="bicor",
                                  maxPOutliers = 0.1,
                                  TOMType = "unsigned", 
                                  minModuleSize = 30,
                                  minCoreKME = 0.3,
                                  minKMEtoStay = 0.1,
                                  deepSplit = 2,
                                  numericLabels = TRUE, 
                                  pamRespectsDendro = TRUE,
                                  saveTOMs = FALSE,
                                  verbose = 1, 
                                  maxBlockSize=3000)
  
  PYR23group2net <- blockwiseModules(PYR23group2, 
                                     power = 6,
                                     networkType = "unsigned",
                                     corType="bicor",
                                     maxPOutliers = 0.1,
                                     TOMType = "unsigned", 
                                     minModuleSize = 30,
                                     minCoreKME = 0.3,
                                     minKMEtoStay = 0.1,
                                     deepSplit = 2,
                                     numericLabels = TRUE, 
                                     pamRespectsDendro = TRUE,
                                     saveTOMs = FALSE,
                                     verbose = 1, 
                                     maxBlockSize=3000)
  
  PYR56group2net <- blockwiseModules(PYR56group2, 
                                     power = 6,
                                     networkType = "unsigned",
                                     corType="bicor",
                                     maxPOutliers = 0.1,
                                     TOMType = "unsigned", 
                                     minModuleSize = 30,
                                     minCoreKME = 0.3,
                                     minKMEtoStay = 0.1,
                                     deepSplit = 2,
                                     numericLabels = TRUE, 
                                     pamRespectsDendro = TRUE,
                                     saveTOMs = FALSE,
                                     verbose = 1, 
                                     maxBlockSize=3000)
  
  
  SSTgroup2net <- blockwiseModules(SSTgroup2, 
                                   power = 6,
                                   networkType = "unsigned",
                                   corType="bicor",
                                   maxPOutliers = 0.1,
                                   TOMType = "unsigned", 
                                   minModuleSize = 30,
                                   minCoreKME = 0.3,
                                   minKMEtoStay = 0.1,
                                   deepSplit = 2,
                                   numericLabels = TRUE, 
                                   pamRespectsDendro = TRUE,
                                   saveTOMs = FALSE,
                                   verbose = 1, 
                                   maxBlockSize=3000)
  
  VIPgroup2net <- blockwiseModules(VIPgroup2, 
                                   power = 6,
                                   networkType = "unsigned",
                                   corType="bicor",
                                   maxPOutliers = 0.1,
                                   TOMType = "unsigned", 
                                   minModuleSize = 30,
                                   minCoreKME = 0.3,
                                   minKMEtoStay = 0.1,
                                   deepSplit = 2,
                                   numericLabels = TRUE, 
                                   pamRespectsDendro = TRUE,
                                   saveTOMs = FALSE,
                                   verbose = 1, 
                                   maxBlockSize=3000)
  
  
  ###Extract module information
  NetworkList <- c("PVgroup1net", "PYR23group1net", "PYR56group1net", "SSTgroup1net", "VIPgroup1net","PVgroup2net", "PYR23group2net", "PYR56group2net", "SSTgroup2net", "VIPgroup2net")
  NetworkDataList <- c("PVgroup1", "PYR23group1", "PYR56group1", "SSTgroup1", "VIPgroup1","PVgroup2", "PYR23group2", "PYR56group2", "SSTgroup2", "VIPgroup2")
  
  for(j in 1:length(NetworkList)){
    currentNet <- get(NetworkList[j])
    prefix <-gsub("net","", NetworkList[j])
    
    currentmoduleLabels <- currentNet$colors
    currentmoduleColors <- labels2colors(currentNet$colors)
    currentMEs <- currentNet$MEs
    currentgeneTree <- currentNet$dendrograms[[1]]
    currentMEs <- moduleEigengenes(get(NetworkDataList[j]), currentmoduleColors)$eigengenes
    currentMEs <- orderMEs(currentMEs)
    currentMEs <- currentMEs[,!colnames(currentMEs)=="MEgrey"]
    names(currentMEs) <- paste0(prefix, names(currentMEs))
    currentGenesClustered <- sum(currentmoduleColors != "grey")
    
    
    assign(paste0(prefix, "moduleLabels"), currentmoduleLabels)
    assign(paste0(prefix, "moduleColors"), currentmoduleColors)
    assign(paste0(prefix, "geneTree"), currentgeneTree)
    assign(paste0(prefix, "MEs"), currentMEs)
    assign(paste0(prefix, "Genes_clustered"), currentGenesClustered)
  }
  
  
  #Save number of modules and genes clustered
  ResultsFile[c(resultsrow:resultsrow+35),1] <- ncol(PVgroup1MEs)
  ResultsFile[c(resultsrow:resultsrow+35),2] <- ncol(PVgroup2MEs)
  ResultsFile[c(resultsrow:resultsrow+35),3] <- ncol(PYR23group1MEs)
  ResultsFile[c(resultsrow:resultsrow+35),4] <- ncol(PYR23group2MEs)
  ResultsFile[c(resultsrow:resultsrow+35),5] <- ncol(PYR56group1MEs)
  ResultsFile[c(resultsrow:resultsrow+35),6] <- ncol(PYR56group2MEs)
  ResultsFile[c(resultsrow:resultsrow+35),7] <- ncol(SSTgroup1MEs)
  ResultsFile[c(resultsrow:resultsrow+35),8] <- ncol(SSTgroup2MEs)
  ResultsFile[c(resultsrow:resultsrow+35),9] <- ncol(VIPgroup1MEs)
  ResultsFile[c(resultsrow:resultsrow+35),10] <- ncol(VIPgroup2MEs)
  
  ResultsFile[c(resultsrow:resultsrow+35),11] <- PVgroup1Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),12] <- PVgroup2Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),13] <- PYR23group1Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),14] <- PYR23group2Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),15] <- PYR56group1Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),16] <- PYR56group2Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),17] <- SSTgroup1Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),18] <- SSTgroup2Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),19] <- VIPgroup1Genes_clustered
  ResultsFile[c(resultsrow:resultsrow+35),20] <- VIPgroup2Genes_clustered
  
  #check rownames line up for meta-networks
  # gsub(rownames(".*-+", "", PVgroup1MEs)) == gsub(".*-+", "", rownames(PYR23group1MEs)) 
  # gsub(rownames(".*-+", "", PVgroup1MEs)) == gsub(".*-+", "", rownames(PYR56group1MEs)) 
  # gsub(rownames(".*-+", "", SSTgroup1MEs)) == gsub(".*-+", "", rownames(VIPgroup1MEs))
  # gsub(rownames(".*-+", "", PVgroup1MEs)) == gsub(".*-+", "", rownames(VIPgroup1MEs))
  # 
  # gsub(rownames(".*-+", "", PVgroup2MEs)) == gsub(".*-+", "", rownames(PYR23group2MEs)) 
  # gsub(rownames(".*-+", "", PVgroup2MEs)) == gsub(".*-+", "", rownames(PYR56group2MEs)) 
  # gsub(rownames(".*-+", "", SSTgroup2MEs)) == gsub(".*-+", "", rownames(VIPgroup2MEs))
  # gsub(rownames(".*-+", "", PVgroup2MEs)) == gsub(".*-+", "", rownames(VIPgroup2MEs))
  
  group1MergedMEs <- cbind(PVgroup1MEs, PYR23group1MEs, PYR56group1MEs, VIPgroup1MEs, SSTgroup1MEs)
  group2MergedMEs <- cbind(PVgroup2MEs, PYR23group2MEs, PYR56group2MEs, VIPgroup2MEs, SSTgroup2MEs)
  
  #indices to exclude intra-cellular connections
  PVgroup1max <- ncol(PVgroup1MEs)
  PYR23group1max <- ncol(PYR23group1MEs) + PVgroup1max
  PYR56group1max <- ncol(PYR56group1MEs) + PYR23group1max
  SSTgroup1max <- ncol(SSTgroup1MEs) + PYR56group1max
  VIPgroup1max <- ncol(VIPgroup1MEs) + SSTgroup1max
  
  PVgroup2max <- ncol(PVgroup2MEs)
  PYR23group2max <- ncol(PYR23group2MEs) + PVgroup2max
  PYR56group2max <- ncol(PYR56group2MEs) + PYR23group2max
  SSTgroup2max <- ncol(SSTgroup2MEs) + PYR56group2max
  VIPgroup2max <- ncol(VIPgroup2MEs) + SSTgroup2max
  
  for(j in 1:length(densitylist)){
    
    threshold <- densitylist[j]
    #control network analysis
    group1MEAdjdata <- corr.test(group1MergedMEs, method = "spearman", adjust = "none")
    
    group1MEAdj <- group1MEAdjdata$r
    group1MEAdj[group1MEAdj < 0] <- 0 
    
    group1MEAdjsig <- group1MEAdjdata$p
    
    group1MEAdj[group1MEAdjsig >= threshold] <- 0
    # threshold <- 1 - densitylist[j]
    # #control network analysis
    # 
    # 
    # group1MEAdj <- stats::cor(group1MergedMEs, method="spearman", use="p")
    # group1MEAdj[group1MEAdj < 0] <- 0
    # group1MEAdj[group1MEAdj < quantile(group1MEAdj, threshold)] <- 0
    
    group1MEAdj[1:PVgroup1max, 1:PVgroup1max] <- 0
    group1MEAdj[(PVgroup1max + 1):PYR23group1max,(PVgroup1max + 1):PYR23group1max] <- 0
    group1MEAdj[(PYR23group1max + 1):PYR56group1max,(PYR23group1max + 1):PYR56group1max] <- 0
    group1MEAdj[(PYR56group1max + 1):SSTgroup1max,(PYR56group1max + 1):SSTgroup1max] <- 0
    group1MEAdj[(SSTgroup1max + 1):VIPgroup1max,(SSTgroup1max + 1):VIPgroup1max] <- 0
    group1Graph <- graph_from_adjacency_matrix(group1MEAdj, mode = "undirected", weighted=TRUE, diag=FALSE)
    
    #####THIS IS FIXED NOW - B/C TRUE WAS QUOTED IN THE ABOVE LINE
    group1Edges <- cbind(get.edgelist(group1Graph), E(group1Graph)$weight)
    
    CTlist <- c("PV", "PYR23", "PYR56", "SST", "VIP")
    CTcombos <- combn(CTlist, 2)
    #Group1 first
    for(k in 1:ncol(CTcombos)){
      CT1 <- CTcombos[1,k]
      CT2 <- CTcombos[2,k]
      
      workingEdges <- tryCatch(group1Edges[group1Edges[,1] %in% c(colnames(get(paste0(CT1,"group1MEs"))), colnames(get(paste0(CT2,"group1MEs")))),], error=function(err) NULL)
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      workingEdges <- tryCatch(workingEdges[workingEdges[,2] %in% c(colnames(get(paste0(CT1,"group1MEs"))), colnames(get(paste0(CT2,"group1MEs")))),], error=function(err) NULL)
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      
      if(is.matrix(workingEdges)==TRUE){
        workingEdges <- workingEdges[!(gsub("_ME.*$","",workingEdges[,1]) == gsub("_ME.*$","",workingEdges[,2])),]
      }
      workingEdgesTOTAL <- tryCatch(nrow(workingEdges), error=function(err) NULL)
      if(is.null(workingEdgesTOTAL)){
        workingEdgesTOTAL <- 0
      }
      
      #normalize for # of genes captured in these eigenegene connections - connection number multiplied by (# of genes involved/total clustered) - or mean/weighted mean of genes connected/average module size
      CT1mods <- data.frame("module"=unique(gsub(paste0(CT1,"group1ME"), "", workingEdges[,1])))
      workingCT1mods <- table(get(paste0(CT1,"group1moduleColors")))
      workingCT1conn <- table(gsub(paste0(CT1,"group1ME"), "", workingEdges[,1]))
      CT1mods$Size <- workingCT1mods[match(CT1mods[,1], names(workingCT1mods))]
      CT1mods$Count <- workingCT1conn[match(CT1mods[,1], names(workingCT1conn))]
      
      CT2mods <- data.frame("module"=unique(gsub(paste0(CT2,"group1ME"), "", workingEdges[,2])))
      workingCT2mods <- table(get(paste0(CT2,"group1moduleColors")))
      workingCT2conn <- table(gsub(paste0(CT2,"group1ME"), "", workingEdges[,2]))
      CT2mods$Size <- workingCT2mods[match(CT2mods[,1], names(workingCT2mods))]
      CT2mods$Count <- workingCT2conn[match(CT2mods[,1], names(workingCT2conn))]
      
      workingEdgesNORMfactor <- paste0(sum(CT1mods$Size), "|", sum(CT2mods$Size))
      
      assign(paste0(CT1, "_", CT2, "group1_EdgesTOTAL"), workingEdgesTOTAL)
      assign(paste0(CT1, "_", CT2, "group1_NORMfactor"), workingEdgesNORMfactor)
    }
    
    ResultsFile[resultsrow,21] <- PV_PYR23group1_EdgesTOTAL
    ResultsFile[resultsrow,22] <- PV_PYR56group1_EdgesTOTAL
    ResultsFile[resultsrow,23] <- PV_SSTgroup1_EdgesTOTAL
    ResultsFile[resultsrow,24] <- PV_VIPgroup1_EdgesTOTAL
    ResultsFile[resultsrow,25] <- PYR23_PYR56group1_EdgesTOTAL
    ResultsFile[resultsrow,26] <- PYR23_SSTgroup1_EdgesTOTAL
    ResultsFile[resultsrow,27] <- PYR23_VIPgroup1_EdgesTOTAL
    ResultsFile[resultsrow,28] <- PYR56_SSTgroup1_EdgesTOTAL
    ResultsFile[resultsrow,29] <- PYR56_VIPgroup1_EdgesTOTAL
    ResultsFile[resultsrow,30] <- SST_VIPgroup1_EdgesTOTAL
    
    ResultsFile[resultsrow,41] <- PV_PYR23group1_NORMfactor
    ResultsFile[resultsrow,42] <- PV_PYR56group1_NORMfactor
    ResultsFile[resultsrow,43] <- PV_SSTgroup1_NORMfactor
    ResultsFile[resultsrow,44] <- PV_VIPgroup1_NORMfactor
    ResultsFile[resultsrow,45] <- PYR23_PYR56group1_NORMfactor
    ResultsFile[resultsrow,46] <- PYR23_SSTgroup1_NORMfactor
    ResultsFile[resultsrow,47] <- PYR23_VIPgroup1_NORMfactor
    ResultsFile[resultsrow,48] <- PYR56_SSTgroup1_NORMfactor
    ResultsFile[resultsrow,49] <- PYR56_VIPgroup1_NORMfactor
    ResultsFile[resultsrow,50] <- SST_VIPgroup1_NORMfactor
    
    #Group2 network
    
    group2MEAdjdata <- corr.test(group2MergedMEs, method = "spearman", adjust = "none")
    
    group2MEAdj <- group2MEAdjdata$r
    group2MEAdj[group2MEAdj < 0] <- 0 
    
    group2MEAdjsig <- group2MEAdjdata$p
    
    group2MEAdj[group2MEAdjsig >= threshold] <- 0
    # group2MEAdj <- stats::cor(group2MergedMEs, method="spearman", use="p")
    # group2MEAdj[group2MEAdj < 0] <- 0
    # group2MEAdj[group2MEAdj < quantile(group2MEAdj, threshold)] <- 0
    # 
    group2MEAdj[1:PVgroup2max, 1:PVgroup2max] <- 0
    group2MEAdj[(PVgroup2max + 1):PYR23group2max,(PVgroup2max + 1):PYR23group2max] <- 0
    group2MEAdj[(PYR23group2max + 1):PYR56group2max,(PYR23group2max + 1):PYR56group2max] <- 0
    group2MEAdj[(PYR56group2max + 1):SSTgroup2max,(PYR56group2max + 1):SSTgroup2max] <- 0
    group2MEAdj[(SSTgroup2max + 1):VIPgroup2max,(SSTgroup2max + 1):VIPgroup2max] <- 0
    group2Graph <- graph_from_adjacency_matrix(group2MEAdj, mode = "undirected", weighted=TRUE, diag=FALSE)
    
    #####THIS IS FIXED NOW - B/C TRUE WAS QUOTED IN THE ABOVE LINE
    group2Edges <- cbind(get.edgelist(group2Graph), E(group2Graph)$weight)
    
    CTlist <- c("PV", "PYR23", "PYR56", "SST", "VIP")
    CTcombos <- combn(CTlist, 2)
    #group2 first
    for(k in 1:ncol(CTcombos)){
      CT1 <- CTcombos[1,k]
      CT2 <- CTcombos[2,k]
      
      workingEdges <- tryCatch(group2Edges[group2Edges[,1] %in% c(colnames(get(paste0(CT1,"group2MEs"))), colnames(get(paste0(CT2,"group2MEs")))),],error=function(err) NULL)
      
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      if(is.matrix(workingEdges)==TRUE){
        workingEdges <- workingEdges[!(gsub("_ME.*$","",workingEdges[,1]) == gsub("_ME.*$","",workingEdges[,2])),]
      }
      
      workingEdges <- tryCatch(workingEdges[workingEdges[,2] %in% c(colnames(get(paste0(CT1,"group2MEs"))), colnames(get(paste0(CT2,"group2MEs")))),], error=function(err) NULL)
      
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      
      workingEdgesTOTAL <- tryCatch(nrow(workingEdges), error=function(err) NULL)
      if(is.null(workingEdgesTOTAL)){
        workingEdgesTOTAL <- 0
      }
      
      workingEdgesNORMfactor <- NA
      
      if(workingEdgesTOTAL != 0){
        #normalize for # of genes captured in these eigenegene connections - connection number multiplied by (# of genes involved/total clustered) - or mean/weighted mean of genes connected/average module size
        CT1mods <- data.frame("module"=unique(gsub(paste0(CT1,"group2ME"), "", workingEdges[,1])))
        workingCT1mods <- table(get(paste0(CT1,"group2moduleColors")))
        workingCT1conn <- table(gsub(paste0(CT1,"group2ME"), "", workingEdges[,1]))
        CT1mods$Size <- workingCT1mods[match(CT1mods[,1], names(workingCT1mods))]
        CT1mods$Count <- workingCT1conn[match(CT1mods[,1], names(workingCT1conn))]
        
        CT2mods <- data.frame("module"=unique(gsub(paste0(CT2,"group2ME"), "", workingEdges[,2])))
        workingCT2mods <- table(get(paste0(CT2,"group2moduleColors")))
        workingCT2conn <- table(gsub(paste0(CT2,"group2ME"), "", workingEdges[,2]))
        CT2mods$Size <- workingCT2mods[match(CT2mods[,1], names(workingCT2mods))]
        CT2mods$Count <- workingCT2conn[match(CT2mods[,1], names(workingCT2conn))]
        
        workingEdgesNORMfactor <- paste0(sum(CT1mods$Size), "|", sum(CT2mods$Size))
      }
      assign(paste0(CT1, "_", CT2, "group2_EdgesTOTAL"), workingEdgesTOTAL)
      assign(paste0(CT1, "_", CT2, "group2_NORMfactor"), workingEdgesNORMfactor)
    }
    
    ResultsFile[resultsrow,31] <- PV_PYR23group2_EdgesTOTAL
    ResultsFile[resultsrow,32] <- PV_PYR56group2_EdgesTOTAL
    ResultsFile[resultsrow,33] <- PV_SSTgroup2_EdgesTOTAL
    ResultsFile[resultsrow,34] <- PV_VIPgroup2_EdgesTOTAL
    ResultsFile[resultsrow,35] <- PYR23_PYR56group2_EdgesTOTAL
    ResultsFile[resultsrow,36] <- PYR23_SSTgroup2_EdgesTOTAL
    ResultsFile[resultsrow,37] <- PYR23_VIPgroup2_EdgesTOTAL
    ResultsFile[resultsrow,38] <- PYR56_SSTgroup2_EdgesTOTAL
    ResultsFile[resultsrow,39] <- PYR56_VIPgroup2_EdgesTOTAL
    ResultsFile[resultsrow,40] <- SST_VIPgroup2_EdgesTOTAL
    
    ResultsFile[resultsrow,51] <- PV_PYR23group2_NORMfactor
    ResultsFile[resultsrow,52] <- PV_PYR56group2_NORMfactor
    ResultsFile[resultsrow,53] <- PV_SSTgroup2_NORMfactor
    ResultsFile[resultsrow,54] <- PV_VIPgroup2_NORMfactor
    ResultsFile[resultsrow,55] <- PYR23_PYR56group2_NORMfactor
    ResultsFile[resultsrow,56] <- PYR23_SSTgroup2_NORMfactor
    ResultsFile[resultsrow,57] <- PYR23_VIPgroup2_NORMfactor
    ResultsFile[resultsrow,58] <- PYR56_SSTgroup2_NORMfactor
    ResultsFile[resultsrow,59] <- PYR56_VIPgroup2_NORMfactor
    ResultsFile[resultsrow,60] <- SST_VIPgroup2_NORMfactor
    
    ResultsFile[resultsrow,61] <- densitylist[j]
    
    resultsrow <- resultsrow + 1
  }
  
  write.csv(ResultsFile, file=paste0(group1, "_vs_", group2,"density thresholded permutations transcriptome-level WGCNA #", min(NumList), " to #", max(NumList), ".csv"))
}

)