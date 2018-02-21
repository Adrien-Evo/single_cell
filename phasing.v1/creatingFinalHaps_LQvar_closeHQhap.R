
LQtaken <- as.numeric(as.vector(unique(hap_df[,1])))

LQvar_closeHQhap_ls <- list()
LQvar_closeHQhap_ls <- lapply(1:length(LQtaken), function(v) { buildingLQ_HQhaps_with_closerHQhaps(v) })

lastLQvar_HQhap <- excludingDuplicates(LQvar_closeHQhap_ls)

buildingLQ_HQhaps_with_closerHQhaps <- function(LQvar) {
  
  finalHapsPerLQvar <- list()
  #LQvar <- 2
  df_perLQ <- hap_df[which(hap_df[,1] == LQtaken[LQvar]),]
  for (nbLines in 1:nrow(df_perLQ)) {
    
    hapsToConnect <- as.numeric(as.vector(unlist(df_perLQ[nbLines, 2:3])))
    hapTags <- as.vector(unlist(df_perLQ[nbLines, 4:5]))
    
    hapTmp <- hapList[[hapsToConnect[1]]]
    fullhapNames <- colnames(hapList[[hapsToConnect[1]]])
    
    if (hapTags[2] == "AD") {
      
      hapTmp <- cbind(hapTmp, hapList[[hapsToConnect[2]]])
      fullhapNames <- c(fullhapNames, colnames(hapList[[hapsToConnect[2]]]))
      
    } else if (hapTags[2] == "BC") {
      
      hapTmp <- cbind(hapTmp, hapList[[hapsToConnect[2]]][c(2,1),])
      fullhapNames <- c(fullhapNames, colnames(hapList[[hapsToConnect[2]]]))
      
    } else if (hapTags[2] == "X") {
      
      hapTmp <- hapTmp
    }
    
    if (hapTags[1] == "AD") {

      varConfigTmp <- c(0,1)
      fullhapNames <- c(LQtaken[LQvar], fullhapNames)

    } else if (hapTags[1] == "BC") {

      varConfigTmp <- c(1,0)
      fullhapNames <- c(LQtaken[LQvar], fullhapNames)
    }

    finalHapsPerLQvar[[nbLines]] <- cbind(varConfigTmp,hapTmp)
    colnames(finalHapsPerLQvar[[nbLines]]) <- fullhapNames
  }
  if (length(finalHapsPerLQvar) == 1) {
    
    toReturn <- finalHapsPerLQvar
  } else {
    toReturn <- mergingHaplotypes(finalHapsPerLQvar)
  }
  return(toReturn[[1]])
}


hapList <- lapply(1:length(hapList), function(x) { creatingHapList(x, lastLQvar_HQhap) })

LQ_HQvarNames <- namesHapListFunction(hapList)
LQ_LQ_var <- setdiff(names(table(subMatrix[,2])), LQ_HQvarNames)

no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores)
cl <- makeCluster(no_cores, type="FORK")
# Initiate cluster
#cl <- makeCluster( mpi.universe.size(), type = 'MPI')
clusterExport(cl, c("hapList", "subMatrix", "subOrder","outputLQPhasing","creatingLQhapsWithTags", "taggingLQvar_HQhap_combination", "LQ_LQ_var"))
finalHapList <- list()

finalHapList <- parLapply(cl=cl, 1:length(LQ_LQ_var), function(varPos) { creatingLQhapsWithTags(LQ_LQ_var[varPos]) })

stopCluster(cl)   


allHapsTagsTogether <- unlist(finalHapList)
allHapsTagsTogether <- allHapsTagsTogether[which(allHapsTagsTogether != 0)]
hap_df <- data.frame(matrix(unlist(strsplit(allHapsTagsTogether, split = "_")), ncol = 5, byrow = T))
#HQlinks_support <- lapply(1:(length(hapList)-1), function(l) { hap_df[which(hap_df[,2] == l & hap_df[,3] == l+1),] })
connectionTag <- buildingHQ_connections(length(hapList),outF=HQ_HQ_connectionLog, t=1)

hapToConnect <- sort(as.numeric(union(hap_df[,2], hap_df[,3])[(union(hap_df[,2], hap_df[,3])) != "X"]))

countH <- 1
tags <- c("AD", "BC", "X")
tagInc <- c()
tagDec <- c()
varToExcludeInc <- list()
varToExcludeDec <- list()

while (countH < length(hapToConnect)) {
  
  print(countH)
  
  sub_df_Inc <- hap_df[which(hap_df[,2] == hapToConnect[countH] & hap_df[,3] == (hapToConnect[countH]+1)),]
  tagInc[countH] <- tags[which(table(sub_df_Inc[,5]) == max(table(sub_df_Inc[,5])))]
  varToExcludeInc[[countH]] <- as.vector(sub_df_Inc[,1][which(sub_df_Inc[,5] != tagInc[countH])])
  
  sub_df_Dec <- hap_df[which(hap_df[,2] == (hapToConnect[countH]+1) & hap_df[,3] == hapToConnect[countH]),]
  tagDec[countH] <- tags[which(table(sub_df_Dec[,5]) == max(table(sub_df_Dec[,5])))]
  varToExcludeDec[[countH]] <- as.vector(sub_df_Dec[,1][which(sub_df_Dec[,5] != tagDec[countH])])
  countH <- countH+1
}

varToExcludeAll <- c(unlist(varToExcludeInc), unlist(varToExcludeDec))
match(varToExcludeAll, hap_df[,1])

hap_df <- hap_df[-match(varToExcludeAll, hap_df[,1]),]

write.table(lastLQvar_HQhap, file=paste("fullHaplotypes_27Aug_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".txt", sep=""), quote = F, row.names = F, col.names = T, sep = "\t")

