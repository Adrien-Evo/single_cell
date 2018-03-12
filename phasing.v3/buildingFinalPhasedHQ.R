require(cluster)

pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/LQvar_HQhap_table_chr10_HQ25_10_5_Mar1"
pathToHQHQ_m <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/HQhap_matrix_chr10_HQ25_10_5_Mar1"
pathToOutputFolder <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/finalHaps_chr10_HQ25_10_5_Mar1"
pathToFinalHapMerging <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/finalHapMerging_chr10_HQ25_10_5_Jan29"
fileNames <- list.files(pathToScratch, pattern="*.txt", full.names=TRUE)

for (file in 7:length(fileNames)) {
  
  #file <- 6
  sufixHQ_HQmatrix <- "HQhapMatrix25_NbCells10_NbLinks5"
  sufixPhasedHap <- "finalHapNeighborHaps_HQ25_NbCells10_NbLinks5"
  idTag <- unlist(strsplit(unlist(strsplit(fileNames[file], split = "/"))[length(unlist(strsplit(fileNames[file], split = "/")))], split="\\."))[1]
  chrTag <- unlist(strsplit(unlist(strsplit(fileNames[file], split = "/"))[length(unlist(strsplit(fileNames[file], split = "/")))], split="\\."))[2]
  cat(paste(idTag, chrTag, sep = "\t"), sep = "\n")
  
  openHQHQ_m <- read.table(paste(pathToHQHQ_m, paste(idTag, chrTag, sufixHQ_HQmatrix, "txt", sep = "."), sep = "/"), header = T)
  namesHQHQ_m <- scan(paste(pathToHQHQ_m, paste(idTag, chrTag, sufixHQ_HQmatrix, "txt", sep = "."), sep = "/"), what = numeric(), nlines = 1)
  colnames(openHQHQ_m) <- namesHQHQ_m
  
  openMatrix <- read.table(fileNames[file], header = F)
  hapList <- creatingHapList(openHQHQ_m)
  #computing distances with gower
  if (ncol(openMatrix) == 2) {
    d <- daisy(matrix(openMatrix[,2:ncol(openMatrix)]), metric = "gower")    
  } else {
    d <- daisy(openMatrix[,2:ncol(openMatrix)], metric = "gower")    
  }

  #replacing the NA entries in the distances by 0.5
  d[is.na(d)] <- 0.5
  #identifying clusters
  clusters <- hclust(d)
  
  countClusters <- 2
  while (countClusters < 15) {
      
    #retrieving indexes of the rows belonging to each cluster
    indx <- cutree(clusters, countClusters)
    #retrieving clusters' ids and amount of var supporting each cluster
    t_indx <- table(indx)
    #retreving the cluster ids supported by the largest amount of variants
    cluster_id <- as.numeric(sort(names(t_indx)[order(t_indx, decreasing = T)[1:2]]))
    
    if (t_indx[cluster_id[1]]/t_indx[cluster_id[2]] < 2) {
      #retrieving the largest configuration tag for each HQ haplotype within a cluster
      l_config_one <- unlist(lapply(2:ncol(openMatrix), FUN = hapConf, cluster=cluster_id[1]))
      l_config_two <- unlist(lapply(2:ncol(openMatrix), FUN = hapConf, cluster=cluster_id[2]))
      
      if (length(l_config_one) > 1) {
        
        config_haps_one <- c(0,slideFunct(as.numeric(names(l_config_one)),1,1))
        config_haps_two <- c(0,slideFunct(as.numeric(names(l_config_two)),1,1))
        if (sum(apply(m <- rbind(config_haps_one, config_haps_two), 2, function(x) { !is.element(x[1], x[2])})) == 0) {
          cat(paste0("Number of clusters found: ", countClusters, "."), sep = "\n")
          break;
        } 
        
      } else if (length(names(l_config_one)) == 1) {
        if (as.numeric(names(l_config_one)) == -(as.numeric(names(l_config_two)))) {
          cat(paste0("Number of clusters found: ", countClusters, "."), sep = "\n")
          break;
        }
      }
    }
    countClusters <- countClusters+1 
  }
  
  if (length(l_config_one) > 1) { # >1 HQ haplotype
    
    config_haps_one <- c(0,slideFunct(as.numeric(names(l_config_one)),1,1))
    config_haps_two <- c(0,slideFunct(as.numeric(names(l_config_two)),1,1))
    
    if (sum(apply(m <- rbind(config_haps_one, config_haps_two), 2, function(x) { !is.element(x[1], x[2])})) == 0) { #checking if major clusters do match
  
      if (names(l_config_one)[1] == "1") {
        
        HQhapOne <- 1
        HQhapTwo <- 2
        
      } else if (names(l_config_one)[1] == "-1") {
        
        HQhapOne <- 2
        HQhapTwo <- 1
      }
      
    
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
      write.table(o_phasedHQhaplotypes[c(which(o_phasedHQhaplotypes[,1] == 0),which(o_phasedHQhaplotypes[,1] == 1)),], file=paste0(pathToOutputFolder,"/", paste(idTag, chrTag, sufixPhasedHap,"txt", sep = ".")), quote = F, row.names = F, col.names = T, sep = "\t")
      
      #merging the LQvar in the third cluster
      #apply(m <- rbind(varLQ[which(!is.na(varLQ))], as.numeric(names(l_config_two[which(!is.na(varLQ))]))), 2, function(x) { is.element(x[1], x[2])})
  
    } else { #major clusters DO NOT match
      
      cat("ERROR: The detected major clusters do not match." , sep = "\n")
    }    
  } else { #the l_config is of size one: THERE IS ONLY ONE HQ hap OR INFO FOR ONLY ONE.
      
    if (as.numeric(names(l_config_one)[1]) == -(as.numeric(names(l_config_two)[1]))) {
      
      if (names(l_config_one)[1] == "1") {
        
        HQhapOne <- 1
        HQhapTwo <- 2
        
      } else if (names(l_config_one)[1] == "-1") {
        
        HQhapOne <- 2
        HQhapTwo <- 1
      }
      
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
      write.table(o_phasedHQhaplotypes[c(which(o_phasedHQhaplotypes[,1] == 0),which(o_phasedHQhaplotypes[,1] == 1)),], file=paste0(pathToOutputFolder,"/", paste(idTag, chrTag, sufixPhasedHap,"txt", sep = ".")), quote = F, row.names = F, col.names = T, sep = "\t")    
      
    } else { #major clusters DO NOT match
      
      cat("ERROR: The detected major clusters do not match.", sep = "\n")
      
    }
      
  }
    
  
  oldPhasing <- read.table(paste0(pathToFinalHapMerging, "/", paste(idTag, chrTag, sufixPhasedHap,"txt", sep = ".")), header = T)   
  colnames(oldPhasing) <- scan(paste0(pathToFinalHapMerging, "/", paste(idTag, chrTag, sufixPhasedHap,"txt", sep = ".")), what = numeric(), nlines = 1)
  cat(paste0("Number of phased sites OLD: ", ncol(oldPhasing), " ."), sep = "\n")  
  cat(paste0("Number of phased sites CLUSTERING: ", ncol(o_phasedHQhaplotypes), " ."), sep = "\n")
  
  subset_oldPhasing <- oldPhasing[ ,intersect(colnames(o_phasedHQhaplotypes),colnames(oldPhasing))]
  subset_o_phasedHQhaplotypes <- o_phasedHQhaplotypes[ ,intersect(colnames(o_phasedHQhaplotypes),colnames(oldPhasing))]
  
  similarityProp <- sum(apply(m <- rbind(subset_o_phasedHQhaplotypes[which(subset_o_phasedHQhaplotypes[,1] == 0),], subset_oldPhasing[which(subset_oldPhasing[,1] == 0),]), 2, function(x) { is.element(x[1], x[2])}))/length(intersect(colnames(o_phasedHQhaplotypes),colnames(oldPhasing)))
  
  if (similarityProp > 0.9) {
    
    cat("The clustering methods is fully compatible with the merging method.", sep = "\n")
    
  } else if (similarityProp < 0.1) {
    
    cat("The clustering methods is fully compatible with the merging method.", sep = "\n")
    
  } else if (similarityProp > 0.1 & similarityProp < 0.9) {
    
    cat("ERROR: The clustering methods is not fully compatible with the merging method.", sep = "\n")
  
  }

}
######################
##
## Functions
##
######################

#creating a hapList from the HQmatrix
creatingHapList <- function(HQm) { 
  
  hp_l <- list()
  
  for (line in 1:nrow(HQm)) {
    
    newLine <- HQm[line,][!is.na(HQm[line,])]
    namesNewLine <-  colnames(HQm)[!is.na(HQm[line,])]
    
    names(newLine) <- namesNewLine
    
    hp_l[[line]] <- newLine
    
  }
  
  return(hp_l)
}
##--------
#----

#flipping HQ haplotypes according to the LQ var information
convertHapList <- function(x) { 
  
  z <- rep(0, length(hapList[[x]])); 
  z[which(hapList[[x]] == 0)] <- 1; 
  names(z) <- names(hapList[[x]]); 
  return(z)
  
}
##--------
#----

# retrieving configuration of HQ haplotypes within the major clusters
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
##-------
#----


slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- prod(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

