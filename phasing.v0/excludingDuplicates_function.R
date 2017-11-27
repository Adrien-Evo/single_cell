#########################
##
## New excludingDuplicates - aug 17, 2017
##
#########################
#it takes a list check the names of all  the elements and see whether there are duplicated names.
#if yes, it means two haplotypes can be merged. 
#it merges the overlapping haplotypes and removed them from the list
#it repeats the process until there are no more duplicated sites in the list. 
excludingDuplicates <- function(l) { #takes as FUNCTIONS: creatingHapList, namesHapListFunction, mergingHaplotypes
  
  l_newHapList_tmp <- l
  
  #l_newHapList_tmp <- uniqueHap
  
  #homogeneizes the list in terms of dataframe and all the HQ haps start with the config 0,1
  l_newHapList_tmp <- lapply(1:length(l_newHapList_tmp), function(x) { creatingHapList(x, l_newHapList_tmp) })
  l_uniqueHapListNames <- lapply(1:length(l_newHapList_tmp), function(x) {colnames(l_newHapList_tmp[[x]])})
  v_uniqueHapListNames <- namesHapListFunction(l_newHapList_tmp)
  
  if (length(l_newHapList_tmp) > 1) {
    while(length(l_newHapList_tmp) > 1 & sum(duplicated(v_uniqueHapListNames)) > 0) { # if there is more than one hap in the list
      
      l_finalHapListOfUniqueHaps <- list()
      newHapCount <- 0    
      hapCount <- 1
      while(length(l_newHapList_tmp) >= 1) { #while the list has more than one element
        
        tmpNames <- l_uniqueHapListNames[[hapCount]]
        l_uniqueHapListNames[[hapCount]] <- 0
        
        indexToCompare <- which(unlist(lapply(1:length(l_uniqueHapListNames), function(y) {sum(is.element(tmpNames, l_uniqueHapListNames[[y]]))})) > 0)
        
        if (length(indexToCompare) > 0) { #first hap of the list does overlap with other haps
          
          hapToExport <- mergingHaplotypes(l_newHapList_tmp[c(hapCount, indexToCompare)])[[1]]
          l_newHapList_tmp <- l_newHapList_tmp[-c(hapCount,indexToCompare)]
          l_uniqueHapListNames <- l_uniqueHapListNames[-c(hapCount,indexToCompare)]
          newHapCount <- newHapCount+1
          
          
        } else { # first hap of the list does not overlap with any other hap
          
          hapToExport <- l_newHapList_tmp[[hapCount]]
          l_newHapList_tmp <- l_newHapList_tmp[-c(hapCount)]
          l_uniqueHapListNames <- l_uniqueHapListNames[-c(hapCount)]
          newHapCount <- newHapCount+1
        }
        l_finalHapListOfUniqueHaps[[newHapCount]] <- hapToExport
        
      } #end of merging
      
      l_newHapList_tmp <- l_finalHapListOfUniqueHaps
      #rm(l_finalHapListOfUniqueHaps)
      l_uniqueHapListNames <- lapply(1:length(l_newHapList_tmp), function(x) {colnames(l_newHapList_tmp[[x]])})
      v_uniqueHapListNames <- namesHapListFunction(l_newHapList_tmp)
    } #no more haps have overlapping sites
    return(l_newHapList_tmp)
    
  } else { #there is one single hap in the list
    
    return(l_newHapList_tmp)
    
  }
  
} 