#creating a list with HQ haplotypes
creatingHQhapList <- function(l_counts_order) { #it takes subMatrix and subOrder and return a matrix with all haplotypes
  
  #l_counts_order <- listHQ[[2]]
  
  for (col in 1:nrow(l_counts_order[[1]])) {
    
    vectorList <- list()
    focalPosFirst <- as.character(l_counts_order[[1]][col,2])
    focalPosSec <- as.character(l_counts_order[[1]][col,3])
    linkSupported <- l_counts_order[[2]][col,][1:2]
    # vectorHap_one <- rep(NA, length(names_subM))
    # vectorHap_two <- rep(NA, length(names_subM))
    
    
    if (sum(is.element(1, linkSupported) & is.element(4, linkSupported))) {
      
      vectorList[[col]] <- rbind(c(0,0),c(1,1))
      
      
    } else if (sum(is.element(2, linkSupported) & is.element(3, linkSupported))) {
      
      vectorList[[col]] <- rbind(c(0,1),c(1,0))
    }
    colnames(vectorList[[col]]) <- c(focalPosFirst,focalPosSec)
    

    if (col == 1) {
      
      newHap <- vectorList[[col]]
      
    } else {
      
    newHap <- merge(newHap, vectorList[[col]], by = focalPosFirst)
    }
    
  }
  newHap <- newHap[,order(as.numeric(colnames(newHap)), decreasing = F)]
  namesNewHap <- colnames(newHap)
  newHap <- as.numeric(newHap[which(newHap[,1] == 0),])
  names(newHap) <- namesNewHap
  
    
  #colnames(m_hap_tmp) <- names_subM
  return(newHap)
  
}
##-------------
#-------


creatingHapListFromLinkCountMatrix <- function(haplotypeMatrix) {
  
  ## converting from matrix to list and removing all the NAs
  list_hap_site_tmp <- list()
  
  for (pos in 1:ncol(haplotypeMatrix)) { #pass from hap matrix to hap list
    
    
    if (length(which(haplotypeMatrix[,pos] == 0)) == 1) { #in case there is only one haplotype covering a site
      
      v.new <- haplotypeMatrix[which(haplotypeMatrix[,pos] == 0),]
      v.new_one <- v.new[which(!is.na(v.new))]
      v.new <- haplotypeMatrix[which(haplotypeMatrix[,pos] == 1),]
      v.new_two <- v.new[which(!is.na(v.new))]
      list_hap_site_tmp[[pos]] <- rbind(v.new_one,v.new_two)
      
    } else { #in case there is >1 haplotype covering a site
      
      v.new <- apply(haplotypeMatrix[which(haplotypeMatrix[,pos] == 0),],2,mean, na.rm=T)
      v.new_one <- v.new[which(!is.na(v.new))]
      v.new <- apply(haplotypeMatrix[which(haplotypeMatrix[,pos] == 1),],2,mean, na.rm=T)
      v.new_two <- v.new[which(!is.na(v.new))]
      list_hap_site_tmp[[pos]] <- rbind(v.new_one,v.new_two)
      
    }
  }
  return(list_hap_site_tmp)
}
##---------------
#------------

#creating a matrix with all overlapping haps 
creatingHapMatrixFromLinkCountMatrix <- function(linkCountM_subset, orderLinkCountM_subset) { #it takes subMatrix and subOrder and return a matrix with all haplotypes
  
  
  
  nrow(linkCountM_subset)
  names_subM <- sort(as.numeric(union(linkCountM_subset[,2],linkCountM_subset[,3])))
  
  for (col in 1:nrow(linkCountM_subset)) {
    
    focalPosFirst <- linkCountM_subset[col,2]
    focalPosSec <- linkCountM_subset[col,3]
    linkSupported <- orderLinkCountM_subset[col,][1:2]
    vectorHap_one <- rep(NA, length(names_subM))
    vectorHap_two <- rep(NA, length(names_subM))
    
    
    if (sum(is.element(1, linkSupported) & is.element(4, linkSupported))) {
      
      vectorHap_one[c(which(names_subM == focalPosFirst),which(names_subM == focalPosSec))] <- c(0,0) 
      vectorHap_two[c(which(names_subM == focalPosFirst),which(names_subM == focalPosSec))] <- c(1,1) 
      
    } else {
      
      vectorHap_one[c(which(names_subM == focalPosFirst),which(names_subM == focalPosSec))] <- c(0,1) 
      vectorHap_two[c(which(names_subM == focalPosFirst),which(names_subM == focalPosSec))] <- c(1,0) 
      
    }
    
    if (col == 1) {
      
      m_hap_tmp <- rbind(vectorHap_one,vectorHap_two)
      
    } else {
      
      m_hap_tmp <- rbind(m_hap_tmp,rbind(vectorHap_one,vectorHap_two))
    }
    
  }
  colnames(m_hap_tmp) <- names_subM
  return(m_hap_tmp)
  
}
##------------------
#------------


convertingIntoMatrix <- function(l, l_names) { #l = hapList; l_names=namesHapList
  
  tmp_matrixSortedHap <- matrix(rep(NA, length(l_names)*(2*length(l))), ncol=length(l_names))
  colnames(tmp_matrixSortedHap) <- l_names
  
  i <- 1
  j <- 2
  countHap <- 1
  while (i < nrow(tmp_matrixSortedHap)) {
    
    tmp_matrixSortedHap[i, match(colnames(l[[countHap]]), colnames(tmp_matrixSortedHap))] <- as.vector(unlist(l[[countHap]][1,]))
    tmp_matrixSortedHap[j, match(colnames(l[[countHap]]), colnames(tmp_matrixSortedHap))] <- as.vector(unlist(l[[countHap]][2,]))
    i <- i+2
    j <- j+2
    countHap <- countHap+1
  }
  
  return(tmp_matrixSortedHap)
}
##--------------
#------------
