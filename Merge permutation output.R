options(stringsAsFactors = FALSE)

setwd("C:/Users/dwigh/OneDrive/Desktop/p-val thresh/")

dirlist <- list.dirs()[-1]


#RUN THIS ONCE
# all jobs ran fully
for(d in dirlist){
  setwd(d)
  filelist <- list.files(pattern=".csv")
  workinglength <- 0
  for(i in 1:length(filelist)){
    workingdata <- read.csv(filelist[i])
    check <- sum(complete.cases(workingdata))
    # print(paste(filelist[i], "cases completed:", check))
    workinglength <- workinglength + check
  }
  print(workinglength)
  setwd("..")
}

for(d in dirlist){
  setwd(d)
  filelist <- list.files(pattern=".csv")
  MergedData <- NULL
  for(i in 1:length(filelist)){
    workingdata <- read.csv(filelist[i])
    rowseq <- seq(0, nrow(workingdata), by=36)[-1]
    
    #fill-in module data into missing columns - data ([,c(2:21)]) is in every 36th row and needs to be pasted into the above 35 rows
    for(j in rowseq){
      workingdata[(j-35):(j-1),c(2:21)] <- workingdata[j,2:21]
    }
    workingdata <- workingdata[!(rowSums(is.na(workingdata))==(ncol(workingdata)-1)),]
    MergedData <- rbind(MergedData, workingdata)
  }
  setwd("..")
  write.csv(file=paste0(gsub("\\.\\/", "",d)," merged dthresh permutations.csv"), MergedData)
  
}
write.csv(MergedData, "Merged permutation data.csv")
