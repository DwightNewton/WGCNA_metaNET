#PITT Tetrad moduel numbers and module-trait analysis (differential WGCNA modules)
library(WGCNA)
library(gridExtra)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)
options(scipen=0)
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential")

load("PITT tetrad differential WGCNA input (normalized counts, no trans, no quantile norm).rData")
for(i in list.files(pattern="auto.rData")){load(i)}

#First, get module # and sizes, a distribution (coloured by module size), and total # genes clustered

MElist <- ls(pattern = "MEs")
colourlist <- ls(pattern="moduleColors")

#see what the range of max module sizes are: around 850 so set max to 900
for(i in 1:length(colourlist)){
  test <- get(colourlist[i])
  test <- test[!test == "grey"]
  print(max(table(test)))
}

#Extract module distributions and numbers
holdingHisto <- as.data.frame(matrix(ncol=1, nrow=1))
names(holdingHisto) <- "Var1"
for(i in 1:length(MElist)){
  #15554 genes used as input
  prefix <- gsub("_.*$", "",MElist[i])
  workingMEs <- get(MElist[i])
  workingColours <- get(colourlist[i])
  
  histoData <- as.data.frame(table(workingColours[!workingColours=="grey"]))
  nmodules <- nrow(histoData)
  
  names(histoData) <- c("Var1", prefix)
  holdingHisto <- merge(holdingHisto, histoData, all=TRUE)
  
}

dim(holdingHisto)

histoMelt <- melt(holdingHisto)
histoMelt$CT <- gsub("Control|MDD|Bipolar|SCHIZ","",histoMelt$variable)
histoMelt$CT <- factor(histoMelt$CT, levels = c("PYR23", "PYR56", "SST", "PV", "VIP"))
histoMelt$Dx <- gsub("PV|PYR23|PYR56|SST|VIP","",histoMelt$variable)
histoMelt$Dx <- factor(histoMelt$Dx, levels = c("Control", "MDD", "Bipolar", "SCHIZ"))

#Generate annotations
holdingHisto <- holdingHisto[-nrow(holdingHisto),]
rownames(holdingHisto) <- holdingHisto$Var1
holdingHisto2 <- t(holdingHisto[,-1])
holdingHisto2 <- as.data.frame(holdingHisto2)
holdingHisto2$Size <- rowSums(!is.na(holdingHisto2))
holdingHisto2$GeneNum <- rowSums(holdingHisto2[,c(1:(ncol(holdingHisto2)-1))], na.rm = TRUE)
holdingHisto2$Group <- rownames(holdingHisto2)
holdingHisto2$CT <- gsub("Control|MDD|Bipolar|SCHIZ","",holdingHisto2$Group)
holdingHisto2$CT <- factor(holdingHisto2$CT, levels = c("PYR23", "PYR56", "SST", "PV", "VIP"))
holdingHisto2$Dx <- gsub("PV|PYR23|PYR56|SST|VIP","",holdingHisto2$Group)
holdingHisto2$Dx <- factor(holdingHisto2$Dx, levels = c("Control", "MDD", "Bipolar", "SCHIZ"))
holdingHisto2 <- holdingHisto2[,c("Size", "GeneNum","CT", "Dx")]



#Facet histogram plot by cell-type and Dx
ggplot(histoMelt, aes(x=value, fill=CT))+
  geom_histogram(color="white", binwidth = 125) +
  scale_fill_manual(values=c("#F8766D", "#b01005", "#00BFC4", "#7CAE00", "#C77CFF"))+
  theme_bw()+
  xlim(c(0,2000))+
  geom_text(data=holdingHisto2, aes(label=Size), x=2000, y=15, hjust = 1, vjust = 1)+
  geom_text(data=holdingHisto2, aes(label=GeneNum), x=0, y=15, hjust = 0, vjust = 1)+
  facet_grid(col=vars(CT), row=vars(Dx), scales="fixed")+
  xlab("Module Size") 

#Module eigengene correlations
expdatalist <- ls(pattern="expData")
clinicaldatalist <- ls(pattern="clinicalData")

rownames(PVclinicalData) <- PVclinicalData$X
rownames(PYR23clinicalData) <- PYR23clinicalData$X
rownames(PYR56clinicalData) <- PYR56clinicalData$X
rownames(SSTclinicalData) <- SSTclinicalData$X
rownames(VIPclinicalData) <- VIPclinicalData$X

setwd("Module-trait correlation plots/")
covars <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Allen brain correlations (Jordan covariates)/mapping.human.csv", row.names = 1)
rownames(covars)
covars <- covars[,!names(covars)=="L4"]

for(i in 1:length(expdatalist)){
  #Organize input data (clinical, mod colours, ExpData)
  prefix <- gsub("expData", "", expdatalist[i])
  celltype <- gsub("Control|MDD|Bipolar|SCHIZ", "", prefix)
  workingclinicaldata <- grep(celltype, clinicaldatalist, value = TRUE)
  workingclinicaldata <- get(workingclinicaldata)
  workingexpdata <- get(expdatalist[i])
  workingColours <- get(colourlist[i])
  
  nSamples = nrow(workingexpdata)
  
  #Subset clinicalData to data of interest: Continuous variables: Age, Sizefactor?, PMI, pH, RIN, then append covariates
  workingclinicaldata <- workingclinicaldata[,c("X", "Age", "sizeFactor" , "PMI", "pH", "RIN")]
  workingclinicaldata <- merge(workingclinicaldata, covars, by.x="X", by.y=0, all.x=TRUE)
  rownames(workingclinicaldata) <- workingclinicaldata$X
  
  #subset clinicalData to current group
  workingclinicaldata <- workingclinicaldata[rownames(workingclinicaldata) %in% rownames(workingexpdata),]
  workingclinicaldata <- workingclinicaldata[,-1]
  
  
  # Recalculate MEs with color labels
  ME0 <- moduleEigengenes(workingexpdata, workingColours)$eigengenes
  ME <- orderMEs(ME0)
  moduleTraitCor <- cor(ME, workingclinicaldata)
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  
  #Re-order rows of heatmap by module size (descending order)
  tempNames <- data.frame(temp1 = gsub("ME", "", rownames(moduleTraitCor)))
  tempColours <- data.frame(temp1 = standardColors())
  tempColours$Num <- rownames(tempColours)
  tempFinal <- merge(tempNames, tempColours, by="temp1", all.x=TRUE)
  tempFinal$Num[tempFinal$temp1=="grey"] <- 0  
  tempFinal$temp1 <- paste0("ME", tempFinal$temp1)
  tempFinal$Num <- as.numeric(tempFinal$Num)
  
  moduleTraitCor <- moduleTraitCor[tempFinal$temp1[order(tempFinal$Num)],]
  moduleTraitPvalue <- moduleTraitPvalue[tempFinal$temp1[order(tempFinal$Num)],]
  
  #Remove grey module - correlations are not useful
  moduleTraitCor <- moduleTraitCor[-1,]
  moduleTraitPvalue <- moduleTraitPvalue[-1,]
  
  #Only colour significant correlations (~r<0.51 for n=19)
  moduleTraitCor[abs(moduleTraitCor) < 0.51] <- 0
  
  pdf(paste0(paste0(prefix, " Module-trait relationships.pdf")))
  # Will display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(workingclinicaldata),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste0(prefix, " Module-trait relationships"))
  dev.off()
  
  
  png(paste0(paste0(prefix, " Module-trait relationships.png")), width = 800, height=800)
  # Will display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(workingclinicaldata),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1,
                 zlim = c(-1,1),
                 main = paste0(prefix, " Module-trait relationships"))
  dev.off()
}

setwd("Module-trait correlation plots/Age correlations/")
#Age only plots
for(i in 1:length(expdatalist)){
  #Organize input data (clinical, mod colours, ExpData)
  prefix <- gsub("expData", "", expdatalist[i])
  celltype <- gsub("Control|MDD|Bipolar|SCHIZ", "", prefix)
  workingclinicaldata <- grep(celltype, clinicaldatalist, value = TRUE)
  workingclinicaldata <- get(workingclinicaldata)
  workingexpdata <- get(expdatalist[i])
  workingColours <- get(colourlist[i])
  
  nSamples = nrow(workingexpdata)
  
  #Subset clinicalData to data of interest: Continuous variables: Age, Sizefactor?, PMI, pH, RIN, then append covariates
  workingclinicaldata <- workingclinicaldata[,c("X", "Age")]
  rownames(workingclinicaldata) <- workingclinicaldata$X
  
  #subset clinicalData to current group
  workingclinicaldata <- workingclinicaldata[rownames(workingclinicaldata) %in% rownames(workingexpdata),]
  
  
  # Recalculate MEs with color labels
  ME0 <- moduleEigengenes(workingexpdata, workingColours)$eigengenes
  ME <- orderMEs(ME0)
  moduleTraitCor <- cor(ME, workingclinicaldata)
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  
  #Re-order rows of heatmap by module size (descending order)
  tempNames <- data.frame(temp1 = gsub("ME", "", rownames(moduleTraitCor)))
  tempColours <- data.frame(temp1 = standardColors())
  tempColours$Num <- rownames(tempColours)
  tempFinal <- merge(tempNames, tempColours, by="temp1", all.x=TRUE)
  tempFinal$Num[tempFinal$temp1=="grey"] <- 0  
  tempFinal$temp1 <- paste0("ME", tempFinal$temp1)
  tempFinal$Num <- as.numeric(tempFinal$Num)
  
  moduleTraitCor <- moduleTraitCor[tempFinal$temp1[order(tempFinal$Num)],]
  moduleTraitPvalue <- moduleTraitPvalue[tempFinal$temp1[order(tempFinal$Num)],]
  
  #Remove grey module - correlations are not useful
  moduleTraitCor <- moduleTraitCor[-1,]
  moduleTraitPvalue <- moduleTraitPvalue[-1,]
  
  #Only colour significant correlations (~r<0.51 for n=19)
  moduleTraitCor[abs(moduleTraitCor) < 0.51] <- 0
  
  pdf(paste0(paste0(prefix, " Module-Age heatmap.pdf")))
  # Will display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(workingclinicaldata),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste0(prefix, " Module-trait relationships"))
  dev.off()
  
  
  png(paste0(paste0(prefix, " Module-Age heatmap.png")), width = 300, height=600)
  # Will display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(workingclinicaldata),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols = rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1,
                 zlim = c(-1,1),
                 main = paste0(prefix, " Module-trait relationships"))
  dev.off()
}




#Generate plots with age, for modules that are significantly correlated
setwd("Age correlations/")
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential/Module-trait correlation plots/Age correlations")
#generate gene-significance X module-membership plots for these modules
for(i in 1:length(expdatalist)){
  #Organize input data (clinical, mod colours, ExpData)
  prefix <- gsub("expData", "", expdatalist[i])
  celltype <- gsub("Control|MDD|Bipolar|SCHIZ", "", prefix)
  workingclinicaldata <- grep(celltype, clinicaldatalist, value = TRUE)
  workingclinicaldata <- get(workingclinicaldata)
  workingexpdata <- get(expdatalist[i])
  workingColours <- get(colourlist[i])

  nSamples = nrow(workingexpdata)

  #Subset clinicalData to data of interest: Continuous variables: Age, Sizefactor?, PMI, pH, RIN, then append covariates
  workingclinicaldata <- as.data.frame(workingclinicaldata[,c("X", "Age")])
  rownames(workingclinicaldata) <- workingclinicaldata$X

  #subset clinicalData to current group
  workingclinicaldata <- workingclinicaldata[rownames(workingclinicaldata) %in% rownames(workingexpdata),]

  # Recalculate MEs with color labels
  ME0 <- moduleEigengenes(workingexpdata, workingColours)$eigengenes
  ME <- orderMEs(ME0)
  ME <- ME[,!names(ME)=="MEgrey"]
  #Gene significance and module membership generation
  AgeDF = as.data.frame(workingclinicaldata$Age);
  names(AgeDF) = "Age"
  # names (colors) of the modules
  modNames = substring(names(ME), 3)
  geneModuleMembership = as.data.frame(cor(workingexpdata, ME, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) <- paste("MM", modNames, sep="")
  names(MMPvalue) <- paste("p.MM", modNames, sep="")
  
  geneTraitSignificance <- as.data.frame(cor(workingexpdata, AgeDF, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  names(geneTraitSignificance) <- paste("GS.", names(AgeDF), sep="")
  names(GSPvalue) <- paste("p.GS.", names(AgeDF), sep="")
  
  for(j in 1:ncol(ME)){
    if(cor.test(ME[,j], AgeDF$Age)$p.value < 0.05){
      
      
      # Define variable weight containing the weight column of datTrait
      
      module <- gsub("ME","", names(ME)[j])
      
      column <- match(module, modNames)
      moduleGenes <- workingColours==module
      
      
      
      png(paste0(prefix, " ", names(ME)[j], " Age correlation.png"), width = 400, height = 400)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for Age",
                         main = paste(prefix,"Module membership vs. gene significance\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      png(paste0(prefix, " ", names(ME)[j], " Age correlation (no labels).png"), width = 400, height = 400)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = "",
                         ylab = "",
                         main = paste(prefix, names(ME)[j], "\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      pdf(paste0(prefix, " ", names(ME)[j], " Age correlation.pdf"))
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for Age",
                         main = paste(prefix,"Module membership vs. gene significance\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      
      #Negative control plots - no consecutive significant modules so just plot the next/previous one
      
      negcon <- j+1
      if(j==ncol(ME)){negcon <- j-1}
      
      module <- gsub("ME","", names(ME)[negcon])
      column <- match(module, modNames)
      moduleGenes <- workingColours==module
      
      
      
      png(paste0(prefix, " ", names(ME)[negcon], " (neg control) Age correlation (no labels).png"), width = 400, height = 400)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = "",
                         ylab = "",
                         main = paste(prefix,"Negative Control\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      png(paste0(prefix, " ", names(ME)[negcon], " (neg control) Age correlation.png"), width = 400, height = 400)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for Age",
                         main = paste(prefix,"Negative Control\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      pdf(paste0(prefix, " ", names(ME)[negcon], " (neg control) Age correlation.pdf"))
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for Age",
                         main = paste(prefix,"Negative Control\n"),
                         abline = TRUE, cex=1.5, lwd=2,
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
    }
  }
}