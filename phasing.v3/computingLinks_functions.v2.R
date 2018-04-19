#######################
##
## New function Aug 18 2017
##
#######################
computingLinks <- function(tmpSite, listOfVar, HQ_genotypeMatrix, link_r, outFileName) { #it takes the geno M filtered for those variants with > minHQCells 
  #takes the GLOBAL variables: minNbCells, chr
  #it does not return anything, it rather creates a file with all links following the qual. conditions
  # link_r can either be LQ_ratio or HQ_ratio
  # tmpSite <- tmpLQSite 
  # listOfVar <- sub_LQvarNames
  # HQ_genotypeMatrix <- tmp_geno
  # outFileName <- LQlinksLog
  
  nbWindows <- 0
  tmpSiteIndex <- which(listOfVar == tmpSite) #it happens that the LQvar is rm from the list because there is < minHQCells covering it. 

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
      
      list_counts_order <- list() #added by April 4
      list_counts_order <- subsettingLinksMatrix(HQlinks_df,link_r) #added by April 4
      return(list_counts_order) #added by April 4
      if (is.null(list_counts_order)) {
        
        cat(paste("Variant: ", tmpSite, "has no meaningful genotype.", sep=" "), file=outFileName, sep="\n", append = T)
        
      }
    } else {
      
      cat(paste("Variant: ", tmpSite, "has no coverage for variants within HQ haplotypes", sep=" "), file=outFileName, sep="\n", append = T)
      #HQlinks_df <- 0 #rmved by April 4
    }
  } else { #there is no LQvar in the list: length(tmpSiteIndex) == 0
    
    cat(paste("Variant:", tmpSite, "does not have enough cells covering.",sep=" "), file=outFileName, sep="\n", append = T)
    #HQlinks_df <- 0 #rmved by April 4
  }
  #return(HQlinks_df)
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
subsettingByQuality_computingLQLinks <- function(tmpLQSite, varNamesInHQhaps, ratio, outLog) {
  #GLOBAL variables: openGenoq, openGeno, varNames, outputLQPhasing, minHQCells
  #calls functions: computingLinks
  
  if (!is.element(tmpLQSite, varNamesInHQhaps)) {
    
    # tmpLQSite <- 24257708
    # varNamesInHQhaps <- varInHQhaps
    # outlog <- LQlinksLog

    LQvarNames <- c(tmpLQSite,varNamesInHQhaps)
    overlappingSites <- intersect(LQvarNames, varNames)
    
    subGenoQ <- openGenoq[match(overlappingSites, varNames),]
    subGeno <- openGeno[match(overlappingSites, varNames),]
    indexHQ <- apply(subGenoQ, 1, function(x) { length(which(x > 20))})
    
    sub_LQvarNames <- LQvarNames[which(indexHQ >= minHQCells)]
    tmp_geno <- subGeno[which(indexHQ >= minHQCells),]
    
    tmpLQvar_HQvar <- computingLinks(tmpLQSite, sub_LQvarNames, tmp_geno, ratio, outLog) #changed April19
    
  } else {
    
    tmpLQvar_HQvar <- 0

  } 
  return(tmpLQvar_HQvar)
}
##---------------------
#--------------