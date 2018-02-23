######################
##Function to merge haplotypes

mergingHaplotypes <- function(newHapList_tmp) { #list_of_haplotypes to merge
  
  newHapList_def <- list()
  indexList_i <- 1
  indexList_j <- 2
  countHaps <- 1
  
  while(indexList_j <= length(newHapList_tmp)) {
    
    z <- is.element(names(newHapList_tmp[[indexList_i]]), names(newHapList_tmp[[indexList_j]])) 
    #print(sum(z))
    if (sum(z) == 0) { #means that uniqueHap[[indexList]] is not contained in the next haplotype
      #this is in the case the first and nd are non-overlapping
      #option A
      #print("option A")
      print(paste("Haplotype:", indexList_i, "is not contained within the haplotype:", indexList_j, sep=" "))
      newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_i]]
      newHapList_def[[countHaps]] <- newHapList_tmp[[indexList_i]]
      if (indexList_j == length(newHapList_tmp)) {
        
        newHapList_def[[countHaps+1]] <- newHapList_tmp[[indexList_j]]
        break;
      }
      indexList_i <- indexList_j
      indexList_j <- indexList_j+1
      countHaps <- countHaps+1
      if (indexList_i == length(newHapList_tmp)) {
        newHapList_def[[countHaps]] <- newHapList_tmp[[indexList_i]]
      }
      
    } else { #the hap is contained in the hap+1
      
      if (sum(z) == length(names(newHapList_tmp[[indexList_i]]))) { #the hap is fully contained in hap+1
        #print("option B")
        newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_j]]
        
      } else if (sum(z) == length(names(newHapList_tmp[[indexList_j]]))) {  #the hap+1 is fully contained in hap
        #print("option C")
        newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_i]]
        
      } else if (sum(z) != length(names(newHapList_tmp[[indexList_i]])) & sum(z) != length(names(newHapList_tmp[[indexList_j]]))) { #hap is partially contained in hap+1
        
        if (sum(z) == 1) { #hap and hap+1 share one single variant
          #if the allele in position X1 in the haplotype indexList_i is the same of the allele in position X1 in the haplotype indexList_j
          if (newHapList_tmp[[indexList_i]][which(z)] == newHapList_tmp[[indexList_j]][match(names(newHapList_tmp[[indexList_i]])[which(z)], names(newHapList_tmp[[indexList_j]]))]) {
            
            newHapList_tmp[[indexList_i]] <- c(newHapList_tmp[[indexList_i]],newHapList_tmp[[indexList_j]][-c(match(names(newHapList_tmp[[indexList_i]])[which(z)], names(newHapList_tmp[[indexList_j]])))])
            #newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_i]][order(as.numeric(names(newHapList_tmp[[indexList_i]])))]
          
          } else { #if the alleles differ between the hap indexList_i and hap indexList_j
            
            altHap_tmp <- rep(0, length(newHapList_tmp[[indexList_j]]))
            altHap_tmp[which(newHapList_tmp[[indexList_j]] == 0)] <- 1 
            names(altHap_tmp) <- names(newHapList_tmp[[indexList_j]])
            newHapList_tmp[[indexList_i]] <- c(newHapList_tmp[[indexList_i]],altHap_tmp[-c(match(names(newHapList_tmp[[indexList_i]])[which(z)], names(altHap_tmp)))])
            
          }
          #newHapList_tmp[[indexList_i]] <- merge(newHapList_tmp[[indexList_i]], subset(newHapList_tmp[[indexList_j]], by=which(z)))
          
        } else if (sum(z) > 1) { #hap and hap+1 share >1 variant
          #print("option E")
          #index_tmp <- 1:length(which(z))
          
          similarity <- apply(rbind(newHapList_tmp[[indexList_i]][which(z)], newHapList_tmp[[indexList_j]][match(names(newHapList_tmp[[indexList_i]])[which(z)], names(newHapList_tmp[[indexList_j]]))]), 2, function(x) {is.element(x[1], x[2]) })
          
          if (sum(similarity) == length(which(z))) {
            
            newHapList_tmp[[indexList_i]] <- c(newHapList_tmp[[indexList_i]],newHapList_tmp[[indexList_j]][-c(match(names(newHapList_tmp[[indexList_i]])[which(z)], names(newHapList_tmp[[indexList_j]])))])
          
          } else if (sum(similarity) == 0) {
            
            altHap_tmp <- rep(0, length(newHapList_tmp[[indexList_j]]))
            altHap_tmp[which(newHapList_tmp[[indexList_j]] == 0)] <- 1 
            names(altHap_tmp) <- names(newHapList_tmp[[indexList_j]])
            newHapList_tmp[[indexList_i]] <- c(newHapList_tmp[[indexList_i]],altHap_tmp[-c(match(names(newHapList_tmp[[indexList_i]])[which(z)], names(altHap_tmp)))])
            
            
          }
          
          #newHapList_tmp[[indexList_i]] <- merge(newHapList_tmp[[indexList_i]][,-c(which(z)[-c(length(which(z)))])], newHapList_tmp[[indexList_j]], by=colnames(newHapList_tmp[[indexList_i]])[which(z)[length(index_tmp)]])
          
        }
        
      }
      
      if (indexList_j == length(newHapList_tmp)) {
        newHapList_def[[countHaps]] <- newHapList_tmp[[indexList_i]]
        break;
      } else {
        indexList_j <- indexList_j+1 ##ADD exception when you have more than one overlap
      }
      
    }
  }
  
  #sorting the variants within haplotypes
  sortedHapList <- lapply(1:length(newHapList_def), function(x) { newHapList_def[[x]][order(as.numeric(names(newHapList_def[[x]])), decreasing = F)] } )
  sortedHapList_v <- lapply(1:length(sortedHapList), function(x) { tmpHap_v <- as.numeric(sortedHapList[[x]])  sortedHapList[[x]][order(as.numeric(names(newHapList_def[[x]])), decreasing = F)] } )
  
  return(sortedHapList)
}
##-------------------


convertingHaplotypesAlwaysToZero <- function(k) {
  
  
  sort_k <- k[order(as.numeric(names(k)), decreasing = F)]
  namesSorted_k <- names(sort_k)
  
  if (sort_k[1] == 0) {
    
    tmpVectorHap <- sort_k
    names(tmpVectorHap) <- namesSorted_k
    
  } else if(k[1] == 1) {
    
    tmpVectorHap <- rep(0, length(sort_k))
    tmpVectorHap[which(tmpVectorHap == 0)] <- 1
    names(tmpVectorHap) <- namesSorted_k

  }
  return(tmpVectorHap)  
  
}






