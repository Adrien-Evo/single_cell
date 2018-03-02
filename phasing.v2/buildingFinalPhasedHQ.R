require(cluster)

<<<<<<< HEAD
pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/LQvar_HQhap_table"
pathToOutputFolder <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/finalHaps"
fileNames <- list.files(pathToScratch, pattern="*.txt", full.names=TRUE)
openMatrix <- read.table(fileNames[2], header = T)
=======
pathToScratch <- getwd()
openMatrix <- read.table(paste(pathToScratch, "/AD385.chr21.LQvar_HQhap_matrix_HQ25_NbCells10_NbLinks5.txt", sep = "/"), header = T)
>>>>>>> 3e6ee59053dde34ace3f9d16bf261487e3537fef
#computing distances with gower 
d <- daisy(openMatrix[,2:ncol(openMatrix)], metric = "gower")
#replacing the NA entries in the distances by 0.5
d[is.na(d)] <- 0.5
#identifying clusters
clusters <- hclust(d)
<<<<<<< HEAD
#plot(clusters)
=======
plot(clusters)
>>>>>>> 3e6ee59053dde34ace3f9d16bf261487e3537fef
#retrieving indexes of the rows belonging to each cluster
indx <- cutree(clusters, 5)
#retrieving clusters' ids and amount of var supporting each cluster
t_indx <- table(indx)
#retreving the cluster ids supported by the largest amount of variants
cluster_id <- as.numeric(sort(names(t_indx)[order(t_indx, decreasing = T)[1:2]]))
#retrieving the largest configuration tag for each HQ haplotype within a cluster
l_config_one <- unlist(lapply(2:ncol(openMatrix), FUN = hapConf, cluster=cluster_id[1]))
l_config_two <- unlist(lapply(2:ncol(openMatrix), FUN = hapConf, cluster=cluster_id[2]))

if (length(l_config_one) > 1) {
  
  config_haps_one <- c(0,slideFunct(as.numeric(names(l_config_one)),1,1))
  config_haps_two <- c(0,slideFunct(as.numeric(names(l_config_two)),1,1))
  
  if (names(l_config_one)[1] == "1") {
    
    HQhapOne <- 1
    HQhapTwo <- 2
    
  } else if (names(l_config_one)[1] == "-1") {
    
    HQhapOne <- 2
    HQhapTwo <- 1
  }
  
<<<<<<< HEAD
  #swaping HQ haplotypes according to the inferred configuration
  l_tmp_hapList <- lapply(which(config_haps_one == -1), FUN = convertHapList)
  hapList[which(config_haps_one == -1)] <- l_tmp_hapList
  #merging HQhap with unlist
  phasedHapList <- unlist(hapList)
  
  LQvar_in_cis <- rep(0,length(openMatrix[which(indx==HQhapOne), 1]))
  LQvar_in_trans <- rep(1,length(openMatrix[which(indx==HQhapTwo), 1]))
  names(LQvar_in_cis) <- openMatrix[which(indx==HQhapOne), 1]
  names(LQvar_in_trans) <- openMatrix[which(indx==HQhapTwo), 1]
  
  phasedHapList <- c(phasedHapList, LQvar_in_cis, LQvar_in_trans)
  
  altHap <- rep(0, length(phasedHapList))
  altHap[which(phasedHapList == 0)] <- 1
  phasedHQhaplotypes <- rbind(phasedHapList, altHap)
  colnames(phasedHQhaplotypes) <- names(phasedHapList)
  
  
  o_phasedHQhaplotypes <- phasedHQhaplotypes[,order(as.numeric(colnames(phasedHQhaplotypes)))]
  write.table(o_phasedHQhaplotypes, file=paste0(pathToOutputFolder,"/AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
  #merging the LQvar in the third cluster
  apply(m <- rbind(varLQ[which(!is.na(varLQ))], as.numeric(names(l_config_two[which(!is.na(varLQ))]))), 2, function(x) { is.element(x[1], x[2])})
=======
  
>>>>>>> 3e6ee59053dde34ace3f9d16bf261487e3537fef
  
  if (sum(apply(m <- rbind(config_haps_one, config_haps_two), 2, function(x) { is.element(x[1], x[2])})) == length(config_haps_one)) {
    
    #CAT SOMETHING
    #cat("The two major cluster are alternatives.")
    
    if (names(l_config_one)[1] == "1") {
      
      clusterZero <- 1
      clusterOne <- 2
      
    } else if (names(l_config_one)[1] == "-1") {
      
      clusterZero <- 2
      clusterOne <- 1
      
    }
    
    phaseLQ_onCis <- rep(0, length(openMatrix[which(indx==clusterZero), 1]))
    names(phaseLQ_onCis) <- openMatrix[which(indx==clusterZero), 1]
    phaseLQ_onTrans <- rep(1, length(openMatrix[which(indx==clusterOne), 1]))
    names(phaseLQ_onTrans) <- openMatrix[which(indx==clusterOne), 1]
    #swaping HQ haplotypes according to the inferred configuration
    l_tmp_hapList <- lapply(which(config_haps_one == -1), FUN = convertHapList)
    hapList[which(config_haps == -1)] <- l_tmp_hapList
    #merging HQhap with unlist
    phasedHapList <- unlist(hapList)
    phasedHapList <- c(phasedHapList, phaseLQ_onCis, phaseLQ_onTrans)
    
    altHap <- rep(0, length(phasedHapList))
    altHap[which(phasedHapList == 0)] <- 1
    phasedHQhaplotypes <- rbind(phasedHapList, altHap, deparse.level = 0)
    colnames(phasedHQhaplotypes) <- names(phasedHapList)
    phasedHQhaplotypes <- phasedHQhaplotypes[,order(as.numeric(colnames(phasedHQhaplotypes)))]
    
    
    #merging the LQvar in the third cluster
    #apply(m <- rbind(varLQ[which(!is.na(varLQ))], as.numeric(names(l_config_two[which(!is.na(varLQ))]))), 2, function(x) { is.element(x[1], x[2])})
    
    
  } else { #clusting analysis did not find major clusters
    
    
    
  }
  
} else { #the l_config is of size one: THERE IS ONLY ONE HQ hap OR INFO FOR ONLY ONE.
  
  if (length(hapList) > 1) {
    
    
    
  } else {
    
    
  }
  
  
}


openHQ_old <- read.table("C:/Users/ialves/Documents/Isabel_OICR/single_cell/phasing.v2/AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt", header = T)
colnames(openHQ_old) <- scan("C:/Users/ialves/Documents/Isabel_OICR/single_cell/phasing.v2/AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt", what = numeric(), nlines = 1)


<<<<<<< HEAD
sum(apply(m <- rbind(o_phasedHQhaplotypes[1,intersect(colnames(o_phasedHQhaplotypes),colnames(openHQ_old))], openHQ_old[1,intersect(colnames(o_phasedHQhaplotypes),colnames(openHQ_old))]), 2, function(x) { is.element(x[1], x[2])}))/length(intersect(colnames(o_phasedHQhaplotypes),colnames(openHQ_old)))
=======
oldPhasing <- read.table("AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5_old.txt", header = T)
newPhasing <- read.table("AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt", header=T)
names(oldPhasing) <- scan("AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5_old.txt", what = numeric(), nlines = 1)

names(newPhasing) <- scan("AD173.chr21.finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt", what = numeric(), nlines = 1)



>>>>>>> 3e6ee59053dde34ace3f9d16bf261487e3537fef

######################
##
## Functions
##
######################

convertHapList <- function(x) { 
  
  z <- rep(0, length(hapList[[x]])); 
  z[which(hapList[[x]] == 0)] <- 1; 
  names(z) <- names(hapList[[x]]); 
  return(z)
  
}


hapConf <- function(x, cluster) { #GLOBAL variable: openMatrix,indx,cluster_id
  
  
  l_one <- table(openMatrix[which(indx==cluster_id[cluster]), x])

  if (length(l_one) > 1) {
    
    tmp_l_one <- l_one[order(l_one, decreasing = T)[1]]
    prop_l_one <- tmp_l_one/sum(l_one)
    

    } else {
    
    prop_l_one <- l_one/sum(l_one)

  }

  return(prop_l_one)
}

#plot(clusters)

slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- prod(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

