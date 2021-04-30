######## ONLY RUN THIS ONCE - PITT TETRAD
#Test speed first on cluster using openBLAS and blockwise n=2000 (single threaded)
#If fast, will do 100k (allowing for better multiple comparison correction) - was not the case, did 10k
########################################################################
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/WGCNA (differential)/Analysis - Differential/Permutation testing/")
set.seed(12345)
PermutationList <- rep(list(NA), 100000)
for(i in 1:100000){
  PermutationList[[i]] <- sample(c(1:35), size = 35, replace = FALSE)
}
save(PermutationList, file="PermutationMasterList.rData", version = 2)

#Random numbers and segment into 200-length-sized pieces (processing chunks)
setwd("..")
load ("PITT tetrad differential WGCNA input (normalized counts, no trans, no quantile norm).rData")
setwd("Permutation testing")

PVmerged <- rbind(PVControlexpData, PVMDDexpData, PVBipolarexpData, PVSCHIZexpData)
PYR23merged <- rbind(PYR23ControlexpData, PYR23MDDexpData, PYR23BipolarexpData, PYR23SCHIZexpData)
PYR56merged <- rbind(PYR56ControlexpData, PYR56MDDexpData, PYR56BipolarexpData, PYR56SCHIZexpData)
SSTmerged <- rbind(SSTControlexpData, SSTMDDexpData, SSTBipolarexpData, SSTSCHIZexpData)
VIPmerged <- rbind(VIPControlexpData, VIPMDDexpData, VIPBipolarexpData, VIPSCHIZexpData)

ControlIDs <- gsub(".*-+","",rownames(PVControlexpData))
MDDIDs <- gsub(".*-+","",rownames(PVMDDexpData))
BipolarIDs <- gsub(".*-+","",rownames(PVBipolarexpData))
SCHIZIDs <- gsub(".*-+","",rownames(PVSCHIZexpData))

#Arrange rownames to be in the same order - insert missing PV library
#Same dimensions and subjects, just need to reorder
rowsequence <- c(gsub(".*-+","",rownames(PVControlexpData)), gsub(".*-+","",rownames(PVMDDexpData)), gsub(".*-+","",rownames(PVBipolarexpData)), gsub(".*-+","",rownames(PVSCHIZexpData)))

PVmerged <- PVmerged[order(match(gsub(".*-+","",rownames(PVmerged)), rowsequence)),]
PYR23merged <- PYR23merged[order(match(gsub(".*-+","",rownames(PYR23merged)), rowsequence)),]
PYR56merged <- PYR56merged[order(match(gsub(".*-+","",rownames(PYR56merged)), rowsequence)),]
SSTmerged <- SSTmerged[order(match(gsub(".*-+","",rownames(SSTmerged)), rowsequence)),]
VIPmerged <- VIPmerged[order(match(gsub(".*-+","",rownames(VIPmerged)), rowsequence)),]

gsub(".*-+","",rownames(PVmerged)) == gsub(".*-+","",rownames(PYR23merged))
gsub(".*-+","",rownames(PVmerged)) == gsub(".*-+","",rownames(PYR56merged))
gsub(".*-+","",rownames(PYR56merged)) == gsub(".*-+","",rownames(SSTmerged))
gsub(".*-+","",rownames(SSTmerged)) == gsub(".*-+","",rownames(VIPmerged))

save(PVmerged, PYR23merged, PYR56merged, SSTmerged, VIPmerged, file="MergedWGCNAFiles_forpermutation.rData", version = 2)
########################################################################
