#converts all the elements of the list in the same format
#aug 16 added the ability to always produce the same type of list
#where the first hap ALWAYS starts with zero
creatingHapList <- function(listOfHaplotypes) {
  
  # listOfHaplotypes <- l_hap
  # nbHap <- 1
  
  thap_one <- as.numeric(listOfHaplotypes[1,])
  thap_two <- as.numeric(listOfHaplotypes[2,])
  
  if (thap_one[1] == 0 & thap_two[1] == 1) {
    
    tmpM <- rbind(as.numeric(listOfHaplotypes[1,]), as.numeric(listOfHaplotypes[2,]))
    colnames(tmpM) <- colnames(listOfHaplotypes)
    
  } else if (thap_one[1] == 1 & thap_two[1] == 0) {
    
    tmpM <- rbind(as.numeric(listOfHaplotypes[2,]), as.numeric(listOfHaplotypes[1,]))
    colnames(tmpM) <- colnames(listOfHaplotypes)
  }
  
  return(tmpM)
}

###-----------
#----------

###############
##Function to get the names of var in the hap list

namesHapListFunction <- function(n) {
  
  varNamesPhasedFirst <- c()
  for (hap in 1:length(n)) {
    
    varNamesPhasedFirst <- c(varNamesPhasedFirst,colnames(n[[hap]]))
    
  }
  return(varNamesPhasedFirst)
}
###-----------
#----------
