######################
##Function to merge haplotypes

mergingHaplotypes <- function(list_of_haplotypes) {
  
  newHapList_tmp <- list_of_haplotypes
  #newHapList_tmp <- newHap
  #newHapList_tmp <- listFakeHaps
  newHapList_def <- list()
  
  
  indexList_i <- 1
  indexList_j <- 2
  countHaps <- 1
  
  while(indexList_j <= length(newHapList_tmp)) {
    
    z <- is.element(colnames(newHapList_tmp[[indexList_i]]), colnames(newHapList_tmp[[indexList_j]])) 
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
      
      if (sum(z) == length(colnames(newHapList_tmp[[indexList_i]]))) { #the hap is fully contained in hap+1
        #print("option B")
        newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_j]]
        
      } else if (sum(z) == length(colnames(newHapList_tmp[[indexList_j]]))) {  #the hap+1 is fully contained in hap
        #print("option C")
        newHapList_tmp[[indexList_i]] <- newHapList_tmp[[indexList_i]]
        
      } else if (sum(z) != length(colnames(newHapList_tmp[[indexList_i]])) & sum(z) != length(colnames(newHapList_tmp[[indexList_j]]))) { #hap is partially contained in hap+1
        
        if (sum(z) == 1) { #hap and hap+1 share one single variant
          #print("option D")
          newHapList_tmp[[indexList_i]] <- merge(newHapList_tmp[[indexList_i]], subset(newHapList_tmp[[indexList_j]], by=which(z)))
          
        } else if (sum(z) > 1) { #hap and hap+1 share >1 variant
          #print("option E")
          index_tmp <- 1:length(which(z))
          newHapList_tmp[[indexList_i]] <- merge(newHapList_tmp[[indexList_i]][,-c(which(z)[-c(length(which(z)))])], newHapList_tmp[[indexList_j]], by=colnames(newHapList_tmp[[indexList_i]])[which(z)[length(index_tmp)]])
          
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
  sortedHapList <- lapply(1:length(newHapList_def), function(x) { newHapList_def[[x]][,order(as.numeric(colnames(newHapList_def[[x]])), decreasing = F)] } ) 
  return(sortedHapList)
}
##-------------------
