#######################
##
##
#######################

buildingHQ_connections <- function(hqListSize, outF, t) { #GLOBAL variable: 1) hap_df.
  
  # hqListSize <- length(hapList)
  # outF <- HQ_HQ_connectionLog
  # t <- minLQsupport
  
  hapCoveredByLQ <- intersect(1:hqListSize, sort(as.numeric(union(hap_df[,2], hap_df[,3]))[!is.na(as.numeric(union(hap_df[,2], hap_df[,3])))]))
  hapNotCovByLQ <- setdiff(1:hqListSize, sort(as.numeric(union(hap_df[,2], hap_df[,3]))[!is.na(as.numeric(union(hap_df[,2], hap_df[,3])))]))
  if (length(hapNotCovByLQ) > 0) {
    cat(paste("Hap Nb:", hapNotCovByLQ, "excluded." , sep=" "), file = outF, append = T, sep = "\n")
  }
  HQhap_index_x <- 1
  HQhap_index_y <- HQhap_index_x+1
  linksTags <- c("AD", "BC", "X")
  tmpConnTag <- c()
  
  while (HQhap_index_x < length(hapCoveredByLQ)) {
    
    # print(hapCoveredByLQ[HQhap_index_x])
    # print(hapCoveredByLQ[HQhap_index_y])
    
    tbl_x <- hap_df[which(hap_df[,2] == hapCoveredByLQ[HQhap_index_x] & hap_df[,3] == hapCoveredByLQ[HQhap_index_y]),]
    
    if (nrow(tbl_x) > 0) { 
      
      tmp_v <- as.vector(table(hap_df[which(hap_df[,2] == hapCoveredByLQ[HQhap_index_x] & hap_df[,3] == hapCoveredByLQ[HQhap_index_y]),5]))
      
      if (sum(tmp_v) >= t & tmp_v[which(tmp_v == max(tmp_v))] >= t*sum(tmp_v[-which(tmp_v == max(tmp_v))])) {
        
        
        tmpConnTag <- c(tmpConnTag, paste(hapCoveredByLQ[HQhap_index_x],  hapCoveredByLQ[HQhap_index_y], linksTags[which(tmp_v == max(tmp_v))], sep = "_"))
        HQhap_index_x <- HQhap_index_y
        HQhap_index_y <- HQhap_index_x+1
        
      } else {
        
        cat(paste("Connection", hapCoveredByLQ[HQhap_index_x], "|", hapCoveredByLQ[HQhap_index_y],  "ERROR: not enough comparisons.", sep=" "), file = outF, append = T, sep = "\n")
        HQhap_index_y <- HQhap_index_y+1
      }
      if (HQhap_index_y > length(hapCoveredByLQ)) {
        break;
      }
    } else {
      
      cat(paste("Connection", hapCoveredByLQ[HQhap_index_x], "|", hapCoveredByLQ[HQhap_index_y],  "has not LQ supporting.", sep=" "), file = outF, append = T, sep = "\n")
      HQhap_index_y <- HQhap_index_y+1
      
      if (HQhap_index_y > length(hapCoveredByLQ)) {
        break;
      }
      
    }
  } 
  return(tmpConnTag)
}



##########################
#####
#### Establishes the connection between each LQ variants and the closest HQ haplotype
#### Generates a list of the hapList length with the anchored LQ variants
#####
##########################

linkingLQvar_toHQhaplotypes <- function() { #takes two GLOBAL variables: hapList, hap_df
  
  tmp_listHQhap_LQvar <- list()
  length(tmp_listHQhap_LQvar) <- length(hapList)
  sitesLQ <- as.numeric(as.vector(unique(hap_df[,1])))
  
  for (siteNb in 1:length(sitesLQ)) {
    
    #siteNb <- 5356
    tmpHQIndex <- unique(as.numeric(as.vector(unlist(hap_df[which(hap_df[,1] == sitesLQ[siteNb]),2:3]))))
    
    if (sum(is.na(tmpHQIndex)) > 0) { #the LQ linked to only one HQ hap
      
      tmpHapNb <- tmpHQIndex[1]
      
    } else { #the LQ linked to >1 HQ hap
      closerHaps <- order(sapply(tmpHQIndex, function(x) { abs(sitesLQ[siteNb]-as.numeric(colnames(hapList[[x]])[1])) } ))
      tmpHapNb <- tmpHQIndex[closerHaps[1]]
      
    }
    
    
    if (length(hap_df[which(hap_df[,1] == sitesLQ[siteNb] & hap_df[,2] == tmpHapNb),4]) != 0) { #the closest HQ hap is in the 2nd col
      
      tmpLinkTag <- hap_df[which(hap_df[,1] == sitesLQ[siteNb] & hap_df[,2] == tmpHapNb),4]
      
    } else if (length(hap_df[which(hap_df[,1] == sitesLQ[siteNb] & hap_df[,3] == tmpHapNb),4]) != 0) { #the closest HQ hap is in the 3rd col
      
      tmpLinkTagDouble <- hap_df[which(hap_df[,1] == sitesLQ[siteNb] & hap_df[,3] == tmpHapNb),4:5]
      if (as.vector(unlist(tmpLinkTagDouble[1])) ==  as.vector(unlist(tmpLinkTagDouble[2]))) {
        tmpLinkTag <- "AD"
      } else {
        tmpLinkTag <- "BC"
      }
      
      
    } else {
      
      print("Error: closest haplotype not found")
      
    }
    
    if (length(tmp_listHQhap_LQvar[[tmpHapNb]]) == 0) {
      tmp_listHQhap_LQvar[[tmpHapNb]] <- paste(sitesLQ[siteNb], tmpLinkTag, sep="_")
    } else {
      tmp_listHQhap_LQvar[[tmpHapNb]] <- c(unlist(tmp_listHQhap_LQvar[[tmpHapNb]]),paste(sitesLQ[siteNb], tmpLinkTag, sep="_"))
    }
    
  }
  return(tmp_listHQhap_LQvar)
}
####-------------------
##-------------


#function to assemble haplotypes
assemblingCompleteHaplotypes <- function(hapNb) {
  
  #hapNb <- 45
  if (length(listHQhap_LQvar[[hapNb]]) > 0) {
    perHQhapLQvar <- data.frame(matrix(unlist(strsplit(listHQhap_LQvar[[hapNb]], split="_")), ncol=2, byrow = T))
    tmpFullHapNames <- c(colnames(hapList[[hapNb]])[1], as.vector(perHQhapLQvar[,1]))
    tmpZeroOneComb <- rep(0,length(as.vector(perHQhapLQvar[,2])))
    tmpZeroOneComb[which(as.vector(perHQhapLQvar[,2]) == "BC")] <- 1
    tmpFullHaplotype <- c(0, tmpZeroOneComb)
    names(tmpFullHaplotype) <- tmpFullHapNames
    
    tmpFullHaplotype <- c(tmpFullHaplotype,hapList[[hapNb]][which(hapList[[hapNb]][,1] == 0),2:ncol(hapList[[hapNb]])])
    
    tmpOppositeHap <- rep(1, length(tmpFullHaplotype))
    tmpOppositeHap[which(tmpFullHaplotype == 1)] <- 0
    names(tmpOppositeHap) <- names(tmpFullHaplotype)
    return(list(tmpFullHaplotype, tmpOppositeHap))
  } else {
    tmpFullHaplotype <- hapList[[hapNb]][1,]
    tmpOppositeHap <- hapList[[hapNb]][2,]
    return(list(tmpFullHaplotype, tmpOppositeHap))
  }  
}

#################################
##
## END OF FUNCTIONS
##
#################################