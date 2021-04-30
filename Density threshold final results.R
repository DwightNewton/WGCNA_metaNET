

####Empirical p-value calculation 
library(reshape2)
library(ggplot2)
library(multtest)
library(cluster)
library(scatterplot3d)
library(WGCNA)
library(dplyr)
library(igraph)
library(psych)
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential/Permutation testing/p-value threshold approach")

mergeddatalist <- list.files(pattern = "merged")
start <- seq(0.1, 0.02, by=-0.01)
final <- c(start, start/10, start/100, start/1000)

for(i in 1:length(mergeddatalist)){
  workingdata <- read.csv(mergeddatalist[i])
  workinggroup <- gsub(" +.*$", "", mergeddatalist[i])
  
  #norm factor didn't generate, will just used total genes (essentially the same anyways)
  workingdata <- workingdata[,-c(1,2,43:62)]
  
  #Perform stats at each density threshold
  densityres <- as.data.frame(matrix(ncol=81, nrow=36))
  densityres[,1] <- final
  
  colnames(densityres) <- c("density", paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_Control_Modules"), 
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_",workinggroup, "_Modules"),  
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_Modules_pval_increase"),
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_Modules_pval_decrease"),
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_Control_GenesClustered"), 
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_",workinggroup, "_GenesClustered"),  
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_GenesClustered_pval_increase"),
  paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"),"_GenesClustered_pval_decrease"),
  paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])),"_Control_connections"), 
  paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "_",workinggroup,"_connections"), 
  paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "_connections_pval_increase"),
  paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "_connections_pval_decrease"))
  
  
  #get real results - at each density threshold
  setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential")
  filelist <- c(list.files(pattern = "Control"), list.files(pattern = workinggroup))
  for(f in filelist){load(f)}
  setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential/Permutation testing/p-value threshold approach")
  densities <- final
  
  #PICK UP HERE
  realdata <- as.data.frame(matrix(ncol=31, nrow=36))
  names(realdata) <- c("density", paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_Module_Diff"),
                       paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_Genes_Diff"),
                       paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "Raw_conDiff"),
                       paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "Norm_conDiff"))
  realdata$density <- densities
  
  ConModulelist <- ls(pattern = paste0("Control_moduleColors"))
  DisModulelist <- ls(pattern = paste0(workinggroup,"_moduleColors"))
  
  #######################################################
  # Changed my mind, made everything normal again (difference in Dis-Con)
  #######################################################
 
  for(j in 1:length(ConModulelist)){
    wcon <- get(ConModulelist[j])
    wdis <- get(DisModulelist[j])
    
    #subtract 1 to account for grey
    wconmod <- length(unique(wcon))-1
    wcongene <- sum(wcon!="grey")
    
    wdismod <- length(unique(wdis))-1
    wdisgene <- sum(wdis!="grey")
    
    moddiff <-  wdismod - wconmod
    genediff <- wdisgene - wcongene
    
    realdata[,j+1] <- moddiff
    realdata[,j+6] <- genediff
  }  
  
 
  
  ConMEs <- ls(pattern = paste0("Control_MEs"))
  DisMEs <- ls(pattern = paste0(workinggroup,"_MEs"))
  for(j in 1:length(ConMEs)){
    wconME <- get(ConMEs[j])
    wconME <- wconME[,colnames(wconME)!="ME0"]
    colnames(wconME) <- paste0(gsub("_MEs","",ConMEs[j]), "_", colnames(wconME))
    assign(ConMEs[j], wconME)
    
    wdisME <- get(DisMEs[j])
    wdisME <- wdisME[,colnames(wdisME)!="ME0"]
    colnames(wdisME) <- paste0(gsub("_MEs","",DisMEs[j]), "_", colnames(wdisME))
    assign(DisMEs[j], wdisME)
  }
  
  #Align rows - critically important
  ConOrder <- gsub(".*Hu", "", rownames(get(ConMEs[1])))
  DisOrder <- gsub(".*Hu", "", rownames(get(DisMEs[1])))
  
  for(j in 1:length(ConMEs)){
    workingConME <- get(ConMEs[j])
    workingConME <- workingConME[match(ConOrder, gsub(".*Hu", "", rownames(workingConME))), ]
    assign(ConMEs[j], workingConME)
    
    workingDisME <- get(DisMEs[j])
    workingDisME <- workingDisME[match(DisOrder, gsub(".*Hu", "", rownames(workingDisME))), ]
    assign(DisMEs[j], workingDisME)
  }
    
  #combine data
  ConMergedMEs <- cbind(get(ConMEs[1]), get(ConMEs[2]), get(ConMEs[3]), get(ConMEs[4]), get(ConMEs[5]))
  DisMergedMEs <- cbind(get(DisMEs[1]), get(DisMEs[2]), get(DisMEs[3]), get(DisMEs[4]), get(DisMEs[5]))
  
  #indices to exclude intra-cellular connections - same order as always: PV, PYR23, PYR56, SST, VIP
  # PVConmax <- ncol(get(ConMEs[1]))
  # PYR23Conmax <- ncol(get(ConMEs[2])) + ncol(get(ConMEs[1]))
  # PYR56Conmax <- ncol(get(ConMEs[3])) + ncol(get(ConMEs[2]))
  # SSTConmax <- ncol(get(ConMEs[4])) + ncol(get(ConMEs[3]))
  # VIPConmax <- ncol(get(ConMEs[5])) + ncol(get(ConMEs[4]))
  # 
  # PVDismax <- ncol(get(DisMEs[1]))
  # PYR23Dismax <- ncol(get(DisMEs[2])) + ncol(get(DisMEs[1]))
  # PYR56Dismax <- ncol(get(DisMEs[3])) + ncol(get(DisMEs[2]))
  # SSTDismax <- ncol(get(DisMEs[4])) + ncol(get(DisMEs[3]))
  # VIPDismax <- ncol(get(DisMEs[5])) + ncol(get(DisMEs[4]))
  
  
  for(j in 1:length(densities)){
    threshold <- densities[j]
    #control network analysis
    ConMEAdjdata <- corr.test(ConMergedMEs, method = "spearman", adjust = "none")
    ConMEAdj <- ConMEAdj <- stats::cor(ConMergedMEs, method="spearman", use="p")
   
    ConMEAdj <- ConMEAdjdata$r
    ConMEAdj[ConMEAdj < 0] <- 0 
      
    ConMEAdjsig <- ConMEAdjdata$p
    
    ConMEAdj[ConMEAdjsig >= threshold] <- 0
    ConGraph <- graph_from_adjacency_matrix(ConMEAdj, mode = "undirected", weighted=TRUE, diag=FALSE)
    

    ConEdges <- cbind(get.edgelist(ConGraph), E(ConGraph)$weight)
    
    CTlist <- c("PV", "PYR23", "PYR56", "SST", "VIP")
    CTcombos <- combn(CTlist, 2)
    #Con first
    for(k in 1:ncol(CTcombos)){
      CT1 <- CTcombos[1,k]
      CT2 <- CTcombos[2,k]
      
      workingEdges <- tryCatch(ConEdges[ConEdges[,1] %in% c(colnames(get(paste0(CT1,"Control_MEs"))), colnames(get(paste0(CT2,"Control_MEs")))),], error=function(err) NULL)
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      workingEdges <- tryCatch(workingEdges[workingEdges[,2] %in% c(colnames(get(paste0(CT1,"Control_MEs"))), colnames(get(paste0(CT2,"Control_MEs")))),], error=function(err) NULL)
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
      
      assign(paste0(CT1, "_", CT2, "Con_EdgesTOTAL"), workingEdgesTOTAL)
    }
    
    #Dis network
    DisMEAdjdata <- corr.test(DisMergedMEs, method = "spearman", adjust = "none")
    
    DisMEAdj <- DisMEAdjdata$r
    DisMEAdj[DisMEAdj < 0] <- 0 
    
    DisMEAdjsig <- DisMEAdjdata$p
    
    DisMEAdj[DisMEAdjsig >= threshold] <- 0
    DisGraph <- graph_from_adjacency_matrix(DisMEAdj, mode = "undirected", weighted=TRUE, diag=FALSE)
    
    #####THIS IS FIXED NOW - B/C TRUE WAS QUOTED IN THE ABOVE LINE
    DisEdges <- cbind(get.edgelist(DisGraph), E(DisGraph)$weight)
    
    CTlist <- c("PV", "PYR23", "PYR56", "SST", "VIP")
    CTcombos <- combn(CTlist, 2)
    #Dis first
    for(k in 1:ncol(CTcombos)){
      CT1 <- CTcombos[1,k]
      CT2 <- CTcombos[2,k]
      
      workingEdges <- tryCatch(DisEdges[DisEdges[,1] %in% c(colnames(get(paste0(CT1, workinggroup, "_MEs"))), colnames(get(paste0(CT2,workinggroup,"_MEs")))),],error=function(err) NULL)
      if(is.null(workingEdges)==FALSE){
        if(is.matrix(workingEdges)==FALSE){
          workingEdges <- matrix(workingEdges, ncol=3, byrow=TRUE)
        }
      }
      workingEdges <- tryCatch(workingEdges[workingEdges[,2] %in% c(colnames(get(paste0(CT1,workinggroup,"_MEs"))), colnames(get(paste0(CT2,workinggroup,"_MEs")))),], error=function(err) NULL)
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
      
      assign(paste0(CT1, "_", CT2, "Dis_EdgesTOTAL"), workingEdgesTOTAL)
    }
    
    ConEdgeData <- ls(pattern="Con_EdgesTOTAL")
    DisEdgeData <- ls(pattern="Dis_EdgesTOTAL")
    
      #Connections, normalized for # of modules
    for(k in 1:10){
      densityres[j,k+41] <- get(ConEdgeData[k])
      densityres[j,k+51] <- get(DisEdgeData[k])
    }
      
    
    
    realdata[j,12] <-  PV_PYR23Dis_EdgesTOTAL - PV_PYR23Con_EdgesTOTAL
    realdata[j,13] <-  PV_PYR56Dis_EdgesTOTAL - PV_PYR56Con_EdgesTOTAL
    realdata[j,14] <-  PV_SSTDis_EdgesTOTAL - PV_SSTCon_EdgesTOTAL
    realdata[j,15] <-  PV_VIPDis_EdgesTOTAL - PV_VIPCon_EdgesTOTAL
    realdata[j,16] <-  PYR23_PYR56Dis_EdgesTOTAL - PYR23_PYR56Con_EdgesTOTAL
    realdata[j,17] <-  PYR23_SSTDis_EdgesTOTAL - PYR23_SSTCon_EdgesTOTAL
    realdata[j,18] <-  PYR23_VIPDis_EdgesTOTAL - PYR23_VIPCon_EdgesTOTAL
    realdata[j,19] <-  PYR56_SSTDis_EdgesTOTAL - PYR56_SSTCon_EdgesTOTAL
    realdata[j,20] <-  PYR56_VIPDis_EdgesTOTAL - PYR56_VIPCon_EdgesTOTAL
    realdata[j,21] <-  SST_VIPDis_EdgesTOTAL - SST_VIPCon_EdgesTOTAL
    
    
    

    
    #Gene-corrected: measure is: number of connections, weighted by: gene/module ratio in each CT (if a greater proportion of co-expressed genes are captured, then the number is higher)
    CTcombos
    reference <- data.frame("CT"=c("PV", "PYR23", "PYR56", "SST", "VIP"), "index"=c(1,2,3,4,5))
    for(k in 1:ncol(CTcombos)){
      z1 <- reference$index[reference$CT==CTcombos[1,k]]
      z2 <- reference$index[reference$CT==CTcombos[2,k]]
      workingcombo <- paste0(CTcombos[1,k],"_",CTcombos[2,k])
      
      ## Multiplicative norm (i.e. the # of connections/number of possible connections)
      Distempcombodata <- get(paste0(workingcombo, "Dis_EdgesTOTAL")) / (ncol(get(DisMEs[z1])) * ncol(get(DisMEs[z2])))
      Contempcombodata <- get(paste0(workingcombo, "Con_EdgesTOTAL")) / (ncol(get(ConMEs[z1])) * ncol(get(ConMEs[z2])))

      realdata[j,k+21] <- (Distempcombodata - Contempcombodata)
    }

    }
  
  
  #get permuted data
  
  densitylist <- final
  
  head(workingdata)
  sum(is.na(workingdata))
  dim(workingdata)
  workingdata <- workingdata[!is.na(workingdata$density),]
  
  dim(workingdata)
  #Subsetting workingdata is being super weird - will split into different densities and save as a list of variable names - needs to be done manually
  densitydata1<-workingdata[workingdata$density == 0.1,]
  densitydata2<-workingdata[workingdata$density == 0.09,]
  densitydata3<-workingdata[workingdata$density == 0.08,]
  densitydata4<-workingdata[workingdata$density == 0.07,]
  densitydata5<-workingdata[workingdata$density == 0.06,]
  densitydata6<-workingdata[workingdata$density == 0.05,]
  densitydata7<-workingdata[workingdata$density == 0.04,]
  densitydata8<-workingdata[workingdata$density == 0.03,]
  densitydata9<-workingdata[workingdata$density == 0.02,]
  densitydata10<-workingdata[workingdata$density == 0.01,]
  densitydata11<-workingdata[workingdata$density == 0.009,]
  densitydata12<-workingdata[workingdata$density == 0.008,]
  densitydata13<-workingdata[workingdata$density == 0.007,]
  densitydata14<-workingdata[workingdata$density == 0.006,]
  densitydata15<-workingdata[workingdata$density == 0.005,]
  densitydata16<-workingdata[workingdata$density == 0.004,]
  densitydata17<-workingdata[workingdata$density == 0.003,]
  densitydata18<-workingdata[workingdata$density == 0.002,]
  densitydata19<-workingdata[workingdata$density == 0.001,]
  densitydata20<-workingdata[workingdata$density == 0.0009,]
  densitydata21<-workingdata[workingdata$density == 0.0008,]
  densitydata22<-workingdata[workingdata$density == 0.0007,]
  densitydata23<-workingdata[workingdata$density == 0.0006,]
  densitydata24<-workingdata[workingdata$density == 0.0005,]
  densitydata25<-workingdata[workingdata$density == 0.0004,]
  densitydata26<-workingdata[workingdata$density == 0.0003,]
  densitydata27<-workingdata[workingdata$density == 0.0002,]
  densitydata28<-workingdata[workingdata$density == 0.0001,]
  densitydata29<-workingdata[workingdata$density == 0.00009,]
  densitydata30<-workingdata[workingdata$density == 0.00008,]
  densitydata31<-workingdata[workingdata$density == 0.00007,]
  densitydata32<-workingdata[workingdata$density == 0.00006,]
  densitydata33<-workingdata[workingdata$density == 0.00005,]
  densitydata34<-workingdata[workingdata$density == 0.00004,]
  densitydata35<-workingdata[workingdata$density == 0.00003,]
  densitydata36<-workingdata[workingdata$density == 0.00002,]
  
  densitylist <- c("densitydata1","densitydata2","densitydata3","densitydata4","densitydata5","densitydata6","densitydata7","densitydata8","densitydata9","densitydata10","densitydata11","densitydata12","densitydata13","densitydata14","densitydata15","densitydata16","densitydata17","densitydata18","densitydata19","densitydata20","densitydata21","densitydata22","densitydata23","densitydata24","densitydata25","densitydata26","densitydata27","densitydata28","densitydata29","densitydata30","densitydata31","densitydata32","densitydata33","densitydata34","densitydata35","densitydata36")
  
  reference2 <- data.frame("CT"=c("PV", "PYR23", "PYR56", "SST", "VIP"), "index"=c(1,3,5,7,9))
  
  for(j in 1:length(densitylist)){
    
    densitydata <- get(densitylist[j])
    
    
    #Hold perm results here
    permdata <- as.data.frame(matrix(ncol=30, nrow=nrow(densitydata)))
    names(permdata) <- c(paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_perm_Module_Diff"),
                         paste0(c("PV", "PYR23", "PYR56", "SST", "VIP"), "_perm_Genes_Diff"),
                         paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "_permRaw_conDiff"),
                         paste0(gsub("_Group2_sigconn", "", names(workingdata[31:40])), "_permNorm_conDiff"))
    
    names(permdata)[1]
    names(densitydata)[2]
    
    #Permuted module and gene diffs
    for(k in 0:4){
      #modules
      permdata[,k+1] <- densitydata[,2*k+2] - densitydata[,2*k+1]
      #Genes
      permdata[,k+6] <- densitydata[,2*k+12] - densitydata[,2*k+11]
      
    }
    
    #Permuted raw and normalized connections
    for(k in 1:10){
      CT1 <- gsub("_+.*$","",names(permdata)[k+10])
      CT2 <- gsub("_+.*$","", sub("^.*?_+","", names(permdata)[k+10]))
      
      
      # #Raw (keep on)
      permdata[,k+10] <- densitydata[,k+30] - densitydata[,k+20]

      #Gene Normalized
      z1 <- reference2$index[reference$CT==CTcombos[1,k]]
      z2 <- reference2$index[reference$CT==CTcombos[2,k]]

      workingcombo <- paste0(CTcombos[1,k],"_",CTcombos[2,k])
      
      #Mult-norm
      ## Multiplicative norm (i.e. the # of connections/number of possible connections)
      Distempcombodata <- densitydata[,k+30] / (densitydata[,1+z1] * densitydata[,1+z2])

      Contempcombodata <- densitydata[,k+20] / (densitydata[,z1] * densitydata[,z2])

      
      permdata[,k+20] <- (Distempcombodata - Contempcombodata)
      
    }
   
    names(densitydata)[10+z2]
    
    #Now, we have "permdata" with the permuted data in the same format (for a particular density threshold), as the observed values in "realdata"
    #Calculate empirical p values for that density - do one-tailed increase and decrease (pull correct one based on observed direction)
    
    
    #Module number and genes clustered
    for(k in 1:length(ConModulelist)){
      wcon <- get(ConModulelist[k])
      wdis <- get(DisModulelist[k])
      
      #subtract 1 to account for grey
      wconmod <- length(unique(wcon))-1
      wcongene <- sum(wcon!="grey")
      
      wdismod <- length(unique(wdis))-1
      wdisgene <- sum(wdis!="grey")
      
      densityres[,k+1] <- wconmod
      densityres[,k+21] <- wcongene
      densityres[,k+6] <- wdismod
      densityres[,k+26] <- wdisgene
    }  
    
    
    #Module number and genes clustered pvals
    for(k in 1:5){
      #"Increase" (one-tailed positive)
      module_workingSum_pos <- sum(realdata[j,k+1] <= permdata[,k])
      module_pval_pos <- (module_workingSum_pos + 1)/(nrow(permdata)+1)
      
      gene_workingSum_pos <- sum(realdata[j,k+6] <= permdata[,k+5])
      gene_pval_pos <- (gene_workingSum_pos + 1)/(nrow(permdata)+1)
      
      
      #"Decrease" (one-tailed negative)
      module_workingSum_neg <- sum(realdata[j,k+1] >= permdata[,k])
      module_pval_neg <- (module_workingSum_neg + 1)/(nrow(permdata)+1)
      
      gene_workingSum_neg <- sum(realdata[j,k+6] >= permdata[,k+5])
      gene_pval_neg <- (gene_workingSum_neg + 1)/(nrow(permdata)+1)
      
      
      densityres[,k+11] <- module_pval_pos
      densityres[,k+16] <- module_pval_neg
      densityres[,k+31] <- gene_pval_pos
      densityres[,k+36] <- gene_pval_neg
    }
    
   
    #connection pvals
    for(k in 1:10){

      # #Normalized
      #Increase
      connection_workingSum_pos <- sum(realdata[j,k+21] <= permdata[,k+20])
      connection_pval_pos <- (connection_workingSum_pos + 1)/(nrow(permdata)+1)

      #Decrease
      connection_workingSum_neg <- sum(realdata[j,k+21] >= permdata[,k+20])
      connection_pval_neg <- (connection_workingSum_neg + 1)/(nrow(permdata)+1)

      #Assign
      densityres[j,k+61] <- connection_pval_pos
      densityres[j,k+71] <- connection_pval_neg
    }
  }
  
  #####################Density threshold figure outputs
  # scatterplot (observed value, coloured by p-value and direction (Red=up, blue=down, grey/white=n.s.))
  #             -permutation ranges (vertical line (min/max value)) 
  #             -reference lines at 0
  
  matrix(names(realdata),ncol=1)
  matrix(names(densityres),ncol=1)
  
  #density threshold plotting object
  
  #NORMALIZED
  densityplotting <- densityres[,-c(2:41)]
  densityplotting <- merge(densityplotting, realdata[,c(1,22:31)])

  matrix(names(densityplotting), ncol=1)

    #add min/max in permuted data
  for(j in 1:length(densities)){
    tempdata <- get(densitylist[j])
    for(k in 1:10){
      prefix <- gsub("Norm_conDiff","",names(densityplotting[k+41]))
      
      #Gene-normalized
      z1 <- reference2$index[reference$CT==CTcombos[1,k]]
      z2 <- reference2$index[reference$CT==CTcombos[2,k]]

      workingcombo <- paste0(CTcombos[1,k],"_",CTcombos[2,k])
      
      # #mult-norm
      Distempcombodata <- tempdata[,k+30] / (tempdata[,1+z1] * tempdata[,1+z2])
      Contempcombodata <- tempdata[,k+20] / (tempdata[,z1] * tempdata[,z2])

      densityplotting[densityplotting$density==densities[j], paste0(prefix, "_min")] <- quantile((Distempcombodata - Contempcombodata), 0.05)
      densityplotting[densityplotting$density==densities[j], paste0(prefix, "_max")] <- quantile((Distempcombodata - Contempcombodata), 0.95)

    }
  }
  ###############To chop-off useless small pvalue comparisons
  densityplotting <- densityplotting[densityplotting$density >= 0.001,]
  ################33
  
  densityplotting_melt <- melt(densityplotting, id.vars = "density")
  #Fix inconsistency
  densityplotting_melt$variable <- gsub("Norm_", "_Norm_", densityplotting_melt$variable)
  densityplotting_melt$variable <- gsub("Raw_", "_Raw_", densityplotting_melt$variable)
  
  densityplotting_melt$CT1 <- gsub("_+.*$", "", densityplotting_melt$variable)
  densityplotting_melt$CT2 <- gsub("_+.*$", "", sub("^.*?_+", "", densityplotting_melt$variable))
  densityplotting_melt$connection <- paste0(densityplotting_melt$CT1, "_", densityplotting_melt$CT2)
  
  
  #Normalized
  densityplotting_melt <- densityplotting_melt[!grepl("connections$", densityplotting_melt$variable),]
  
  for(j in 1:nrow(densityplotting_melt)){
    densityplotting_melt$variable[j] <- gsub(paste0(densityplotting_melt$connection[j],"_"), "", densityplotting_melt$variable[j] )
  }
  densityplotting_cast <- dcast(densityplotting_melt, density+CT1+CT2+connection~variable)
  
  #Colour coding points by significance
  
  #NORMALIZED
  densityplotting_cast$Direction <- NA
  for(j in 1:nrow(densityplotting_cast)){
    densityplotting_cast$Direction[j] <- 0

    if((densityplotting_cast$Norm_conDiff[j] < 0) & (densityplotting_cast$connections_pval_decrease[j] < 0.05)){
      densityplotting_cast$Direction[j] <- -1
    }

    if((densityplotting_cast$Norm_conDiff[j] < 0) & (densityplotting_cast$connections_pval_decrease[j] > 0.05) & (densityplotting_cast$connections_pval_decrease[j] < 0.1)){
      densityplotting_cast$Direction[j] <- -0.5
    }
    if((densityplotting_cast$Norm_conDiff[j] > 0) & (densityplotting_cast$connections_pval_increase[j] < 0.05)){
      densityplotting_cast$Direction[j] <- 1
    }

    if((densityplotting_cast$Norm_conDiff[j] > 0) & (densityplotting_cast$connections_pval_decrease[j] > 0.05) & (densityplotting_cast$connections_pval_decrease[j] < 0.1)){
        densityplotting_cast$Direction[j] <- 0.5
    }
  }

  densityplotting_cast$Direction <- as.factor(densityplotting_cast$Direction)
  densityplotting_cast$Ordering <- rep(1:19, each=10)
  densityplotting_cast$density <- as.numeric(densityplotting_cast$density)
  
  #facet_grid diagram showing each ct-specific density threshold pval plot (CT's as both rows/columns)
  #NORMALIZED
  outputplots_grid <- ggplot(densityplotting_cast, aes(x=density, y=Norm_conDiff, fill=Direction)) +
    scale_fill_manual(values=c("-1"="dodgerblue3", "-0.5"="#76B2EF", "0"="white", "0.5"="#FF8C8C", "1"="firebrick3"))+
    geom_segment(aes(x=density, y=min, xend=density, yend=max), colour="black")+
    geom_point(shape=21, size=2, colour="black")+
    geom_hline(aes(yintercept=0), size=1.2)+
    geom_vline(aes(xintercept=0.05), colour="grey60", alpha=0.8, linetype="dashed", size=0.8) +
    scale_x_log10() +
    theme_bw()+
    facet_grid(rows=vars(CT1), cols = vars(CT2)) +
    ylab("Difference in Proportion of Connections") +
    xlab("p-value Threshold")

   
  pdf(paste0(workinggroup," (mult-norm) pval threshold pvalue plots - grid.pdf"))
  print(outputplots_grid)
  dev.off()
  
  write.csv(densityres, file=paste0(workinggroup, " (mult-norm, correct gene p) pval-thresh p-values .csv"))
}









