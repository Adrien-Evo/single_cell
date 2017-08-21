#######################
##
## New function Aug 18 2017
##
#######################
computingLinks <- function(tmpSite, listOfVar, HQ_genotypeMatrix, outFileName) { #it takes the geno M filtered for those variants with > minHQCells 
  #takes the GLOBAL variables: minNbCells, chr
  #it does not return anything, it rather creates a file with all links following the qual. conditions
  
  nbWindows <- 0
  #print(tmpSite)
  tmpSiteIndex <- which(listOfVar == tmpSite)
  #listOfVar <- sub_LQvarNames
  #tmpSiteIndex <- which(listOfVar == tmpLQSite)
  #tmpLQSite, sub_LQvarNames, tmp_geno, outputLQPhasing
  
  if (length(tmpSiteIndex) > 0) {
    
    count_var_y <- tmpSiteIndex+1
    tmp_df_links <- list()
    
    while (count_var_y <= nrow(HQ_genotypeMatrix)) {
      
      tmp <- HQ_genotypeMatrix[c(tmpSiteIndex,count_var_y),]
      
      col_i <- c()
      indexToTake <- as.vector(which(apply(tmp, 2, function(x) { col_i <- c(col_i, sum(is.na(x))) }) == 0))
      #print(paste("Nb of overlapping cells:", length(indexToTake)))
      
      if (length(indexToTake) >= minNbCells) {
        
        
        tmp_sub <- tmp[,indexToTake]
        v_countLinks <- apply(vapply(1:ncol(tmp_sub), get_LinkCount, tmp=tmp_sub, FUN.VALUE = rep(0,4)), 1, sum)
        #print(listOfVar[count_var_y])
        
        #cat("Link found", file=outFileName, sep="\n", append = T)
        nbWindows <- nbWindows+1
        #cat(paste("Link nb:", nbWindows, sep=" "), file=outFileName, sep="\n", append=T)
        #cat(paste("Two SNP window:", listOfVar[tmpSiteIndex], "|", listOfVar[count_var_y], sep=" "), file=outFileName, sep="\n", append = T)
        
        tmp_df_links[[nbWindows]] <- c(chr, tmpSite, listOfVar[count_var_y], v_countLinks[1], v_countLinks[2], v_countLinks[3], v_countLinks[4])
        
      } 
      count_var_y <- count_var_y+1
    }
    if (length(tmp_df_links) > 0) {
      
      HQlinks_df <- data.frame(matrix(unlist(tmp_df_links), ncol=7, byrow = T))
      names(HQlinks_df) <- c("chrNb","PosOne", "PosTwo", "0/0", "0/1", "1/0", "1/1")
    } else {
      HQlinks_df <- 0
    }
  } else {
    
    HQlinks_df <- 0
  }
  return(HQlinks_df)
}
##---------
#------

#Counting links
get_LinkCount <- function(x, tmp) {
  
  count_zz <- 0
  count_zo <- 0
  count_oz <- 0
  count_oo <- 0
  
  if (tmp[1,x] == 0 &  tmp[2,x] == 0) {
    count_zz <- 1
  } else if (tmp[1,x] == 0 &  tmp[2,x] == 1) {
    count_zo <- 1
  } else if (tmp[1,x] == 1 &  tmp[2,x] == 0) {
    count_oz <- 1 
  } else if (tmp[1,x] == 1 &  tmp[2,x] == 1) {
    count_oo <- 1
  }
  return(countsTotal=c(count_zz,count_zo,count_oz,count_oo))
  
}

##-----------
#------
#computed the pairwise combinations for a given LQ var and all the var in the HQ haplotypes
subsettingByQuality_computingLQLinks <- function(tmpLQSite, varNamesInHQhaps) {
  #GLOBAL variables: openGenoq, openGeno, varNames, outputLQPhasing, minHQCells
  #calls functions: computingLinks
  
  if (!is.element(tmpLQSite, varNamesInHQhaps)) {
    
    #tmpLQSite <- 24235027
    #varNamesInHQhaps <- varInHQhaps
    
    LQvarNames <- c(tmpLQSite,varNamesInHQhaps)
    overlappingSites <- intersect(LQvarNames, varNames)
    
    subGenoQ <- openGenoq[match(overlappingSites, varNames),]
    subGeno <- openGeno[match(overlappingSites, varNames),]
    indexHQ <- apply(subGenoQ, 1, function(x) { length(which(x > 20))})
    
    sub_LQvarNames <- LQvarNames[which(indexHQ >= minHQCells)]
    tmp_geno <- subGeno[which(indexHQ >= minHQCells),]
    
    tmpLQvar_HQvar <- computingLinks(tmpLQSite, sub_LQvarNames, tmp_geno, outputLQPhasing)
    
    
    
  } else {
    
    tmpLQvar_HQvar <- 0
    
    
  } 
  return(tmpLQvar_HQvar)
}
##---------------------
#--------------