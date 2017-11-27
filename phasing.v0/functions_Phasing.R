###################################
##
## Functions for phasingPipeline_Links_all.R
##
##################################

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
####-----------------
##--------------
#-----



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

#------
#subsets links matrix according to QC we specify
subsettingLinksMatrix <- function(PC_tbl) { #local variable: PC_tbl = matrix with HQ link; global variables: 1)HQ_ratio, 2) minNbLinks; FUNCTIONS: 1) countingHetLinks.
  
  #PC_tbl <- fullHQlinksComb
  #keep pwise Comb with link counts > minNbLinks
  linksM <- lapply(1:nrow(PC_tbl[,4:7]), function(x) { which(as.numeric(as.matrix(PC_tbl[x,4:7])) >= minNbLinks) })
  vvv <- which(unlist(lapply(linksM, FUN=length)) == 2) #vector containing those positions with two most prevalent links
  
  #write.table(kkk[vvv,], file="clean_AD393.chr15.pairwiseComb.phase2.txt", quote = F, row.names = F, col.names = F)
  
  #subsetting the original table with pwise Comb counts
  matrix_counts_links_ltTWO <- PC_tbl[vvv,]
  
  #computing a table with the link count order 2,3,4,1 means that the hights nb of counts if in column 2 and 3
  order_m <- t(apply(matrix_counts_links_ltTWO[,4:7],1, function(o) { order(o, decreasing = T) })) #order of the most prevalent links among sites in vvv
  #getting the row indexes compatible with SNP het status
  tmpOneColHet <- countingHetLinks(order_m)
  tmpTwoColHet <- which(lapply(1:nrow(matrix_counts_links_ltTWO), function(x) { sum(as.numeric(as.matrix(matrix_counts_links_ltTWO[x,order_m[x,1:2]+3]))) > HQ_ratio*sum(as.numeric(as.matrix(matrix_counts_links_ltTWO[x,order_m[x,3:4]+3]))) }) == T)
  colHet <- intersect(tmpOneColHet,tmpTwoColHet)
  
  #building the m_hap matrix, it includes all the haplotypes that and needs to be clean to rm duplicates and merge overlapping haplotypes
  subMatrix <- matrix_counts_links_ltTWO[colHet,]
  #match(sort(subMatrix[1,3:6], decreasing=T), subMatrix[1,3:6])
  subOrder <- order_m[colHet,]
  
  return(list(subMatrix, subOrder))
  
}
##--------------------
#--------------

##--------------------
#--------------
#this function takes the HQ filtered matrix with the link counts order and 
#checks whether the two highest links indicate a hetSNP AS THEY SHOULD
countingHetLinks <- function(linkCountM) { #it takes the order matrix and it returns a vector with row indexes
  
  colHet <- c()
  
  for (col in 1:nrow(linkCountM)) {
    
    index_Supported_Links <- linkCountM[col,1:2]
    #print(index_Supported_Links)
    
    if (is.element(1,index_Supported_Links) & is.element(4,index_Supported_Links)) {
      
#      print("Both first and second pos are heterozygous")
      colHet <- c(colHet,col)
      
    } else if (is.element(2,index_Supported_Links) & is.element(3,index_Supported_Links)) {
      
#      print("Both first and second pos are heterozygous") 
      colHet <- c(colHet,col)
      
    } 
    
  } 
  return(colHet)
}
##------------------
#------------
#converts all the elements of the list in the same format
#aug 16 added the ability to always produce the same type of list
#where the first hap ALWAYS starts with zero
creatingHapList <- function(nbHap, listOfHaplotypes) {
  
  # listOfHaplotypes <- l_hap
  # nbHap <- 1
  
  thap_one <- as.numeric(listOfHaplotypes[[nbHap]][1,])
  thap_two <- as.numeric(listOfHaplotypes[[nbHap]][2,])
  
  if (thap_one[1] == 0 & thap_two[1] == 1) {
    
    tmpM <- rbind(as.numeric(listOfHaplotypes[[nbHap]][1,]), as.numeric(listOfHaplotypes[[nbHap]][2,]))
    colnames(tmpM) <- colnames(listOfHaplotypes[[nbHap]])
    
  } else if (thap_one[1] == 1 & thap_two[1] == 0) {
    
    tmpM <- rbind(as.numeric(listOfHaplotypes[[nbHap]][2,]), as.numeric(listOfHaplotypes[[nbHap]][1,]))
    colnames(tmpM) <- colnames(listOfHaplotypes[[nbHap]])
  }
  
  return(tmpM)
}


##------------------
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

#Counting links
get_LinkCount <- function(x) {
  
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
###---------------
gettingHetSNPdensity <- function(x) { #x-list from tmp; b-bins; ws-window size
  
  nbSNPS_tmp <- c()
  countBins <- 1
  
  while (countBins < length(bins)) {
    
    if (countBins == (length(bins)-1)) {
      
      nbSNPS_tmp[countBins] <- length(which(x >= bins[countBins] & x <= bins[countBins+1]))
      
    } else {
      
      nbSNPS_tmp[countBins] <- length(which(x >= bins[countBins] & x < bins[countBins+1]))
      
    }
    
    if (nbSNPS_tmp[countBins] == 0) {
      
      nbSNPS_tmp[countBins] <- 0
      
    } else {
      
      nbSNPS_tmp[countBins] <- nbSNPS_tmp[countBins]
    }
    
    # print(nbSNPS[countBins])
    # 
    # cat(paste(chrms, bins[countBins],  bins[countBins+1]-1, nbSNPS[countBins], sep="\t"),file=outFileName, sep="\n", append = T)
    countBins <- countBins+1
    
  }
  return(nbSNPS_tmp)
}

######################
##Function to merge haplotypes

mergingHaplotypes <- function(hapList) {
  
newHapList_tmp <- hapList
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


###############
##Function to get the names of var in the hap list

namesHapListFunction <- function(n) {
  
  varNamesPhasedFirst <- c()
  for (hap in 1:length(n)) {
    
    varNamesPhasedFirst <- c(varNamesPhasedFirst,colnames(n[[hap]]))
    
  }
  return(varNamesPhasedFirst)
}
##-------------------


###############
## Function to read link count tables. it
## takes both the HQ and low Q table with different formats

openLinkCountTable <- function(h, chr, fileName, space) { #h = header(T/F), chr = is there chr #? (T/F); fileName = name of the file with Links; space = type of sep 
  #this function requires global variable minNbLinks
  openPC <- read.table(fileName, header = h, sep=space) 
  
  
  if (chr) {
    
    pairWCpos <- openPC[,2:3]
    counts_m <- as.matrix(openPC[,4:7]) 
    kkk <- t(apply(counts_m, 1, function(z) { as.numeric(unlist(strsplit(z[1:4], split="/"))[c(3,6,9,12)]) }))
    #largerThan10 <- apply(kkk, 1, function(x) { which(x>=minNbLinks)} )
    tableC <- cbind(pairWCpos,kkk)
    
  } else {
    
    pairWCpos <- openPC[,1:2]
    counts_m <- as.matrix(openPC[,3:6]) 
    #largerThan10 <- apply(counts_m, 1, function(x) { which(x>=minNbLinks)} )
    tableC <- cbind(pairWCpos,counts_m)
    
  }
  
  return(tableC)
  
}
##-------------------


computePCshownByMoreThanXCells <- function(m) { #m is the matrix with the counts of the links across pairwise comparisons
  #this function requires GLOBAL variable: minNbLinks
  largerThan10 <- apply(m, 1, function(x) { which(x>=minNbLinks)} )
  
  return(largerThan10)
}
##-------------------

##-------------------
#--------------
#computed the pairwise combinations for a given LQ var and all the var in the HQ haplotypes
subsettingByQuality_computingLQLinks <- function(tmpLQSite, varNamesInHQhaps) {
  #GLOBAL variables: openGenoq, openGeno, varNames, outputLQPhasing, minHQCells
  #calls functions: computingLinks
  
  if (!is.element(tmpLQSite, varNamesInHQhaps)) {
    
    # tmpLQSite <- hetSNPsPos[1]
    # varNamesInHQhaps <- varInHQhaps
    
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

#######################
## NEW FUNCTIONS - july 17 2017
##
#######################
#this function takes the subMatrix following QC and generates the possible haplotypes between a provided LQ variant and all the possible
#and previously built HQ haplotypes. It exports a list with the haplotype
#it takes global variables: hapList, subMatrix, LQhapQC
creatingLQhaplotypes <- function(LQVarCoord) {
  
  #pos <- 3
  #LQVarCoord <- "24152993"
  #checks how many pairwise comb in the LQ matrix are related to HQ haplotypes and for how many sites within each HQ hap are covered by the LQ list
  #overlappingHap is a vector of the length of the hapList and the numbers represent the nb of sites with links information in the LQ matrix
  overlappingHap <- unlist(lapply(1:length(hapList), function(x) { sum(is.element(subMatrix[which(subMatrix[,1] == LQVarCoord),2], colnames(hapList[[x]]))) }))
  indexOverlappingHap <-  which(overlappingHap >= 2)
  print(length(indexOverlappingHap))
  
  if (length(indexOverlappingHap) > 0) {
    
    corrPC_one <- c()
    corrPC_two <- c()
    newHap <- list()  
    #hapNb <- 75
    countHapOver <- 1
    for (hapNb in indexOverlappingHap) {
      
      #subsetting the subMatrix per LQ variant
      perHQhap_NewLQPos_m <- subMatrix[intersect(which(subMatrix[,1] == LQVarCoord),which(subMatrix[,2] %in% as.numeric(colnames(hapList[[hapNb]])))),]
      #subsetting the matrix with the links' order
      perHQhap_NewLQPos_order_m <- subOrder[intersect(which(subMatrix[,1] == LQVarCoord),which(subMatrix[,2] %in% as.numeric(colnames(hapList[[hapNb]])))),]
      vectorList <- list()
      
      
      for (col in 1:nrow(perHQhap_NewLQPos_m)) {
        
        linkSupported <- perHQhap_NewLQPos_order_m[col,][1:2]
        
        
        
        if (is.element(1, linkSupported) & is.element(4, linkSupported)) {
          
          vectorList[[col]] <- rbind(c(0,0),c(1,1))
          colnames(vectorList[[col]]) <- c(perHQhap_NewLQPos_m[col,1], perHQhap_NewLQPos_m[col,2])
          
        } else if (is.element(2, linkSupported) & is.element(3, linkSupported)) {
          
          vectorList[[col]] <- rbind(c(0,1),c(1,0))
          colnames(vectorList[[col]]) <- c(perHQhap_NewLQPos_m[col,1], perHQhap_NewLQPos_m[col,2])
          
        }
        
        if (col == 1) {
          
          newHap[[countHapOver]] <- vectorList[[col]]
          
        } else {
          
          newHap[[countHapOver]] <- merge(newHap[[countHapOver]], vectorList[[col]], by = LQVarCoord)
          
          
        }
        
      }
      
      # the first which does this - looking for the row number with a zero in the first overlapping site between the HQ hap and the new LQ hap that is not the LQ SITE
      #the following vector contain the combination of states followed by the first zero in the HQ haplotype
      hapZero_One <- as.vector(hapList[[hapNb]][which(hapList[[hapNb]][,colnames(hapList[[hapNb]])[colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])][1]] == 0 ),colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])])
      #the following vector contain the combination of states followed by the first zero in the LQ haplotype
      hapZero_Two <- as.vector(newHap[[countHapOver]][which(newHap[[countHapOver]][,colnames(newHap[[countHapOver]])[colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])][1]] == 0 ),colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])])
      #example:
      #HQ hap 
      # SNP1 SNP2
      # 0    1
      # 1    0
      #SNP1 and 2 are overlapping between the HQ hap and the new LQ hap
      # hapZero_One contains coordinate c(0, column index(SNP1, SNP2))
      
      #the two vectors below contain the same as before but the combination starting with 1.
      hapOne_One <- as.vector(hapList[[hapNb]][which(hapList[[hapNb]][,colnames(hapList[[hapNb]])[colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])][1]] == 1 ),colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])])
      hapOne_Two <- as.vector(newHap[[countHapOver]][which(newHap[[countHapOver]][,colnames(newHap[[countHapOver]])[colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])][1]] == 1 ),colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])])
      
      #here we compare the HQ hap with the LQ hap for the overlapping sites
      corrPC_one[countHapOver] <- sum(apply(m <- rbind(hapZero_One, hapZero_Two), 2, function(x) { is.element(x[1], x[2])}))
      corrPC_two[countHapOver] <- sum(apply(m <- rbind(hapOne_One, hapOne_Two), 2, function(x) { is.element(x[1], x[2])}))
      countHapOver <- countHapOver+1
      
    } #close for across compatible HQ haps
    
    #corrPC_one/two contain the nb of EQUAL comparisons between the new and the HQ haplotype
    #if the nb of comparisons is >= 5 the LQ is kept and the LQ haploytpe built
    if (LQhapQC == "loose" & (sum(overlappingHap[overlappingHap >= 2]) >= 5)) {
      if ((sum(corrPC_one)/sum(overlappingHap[overlappingHap >= 2]) == 1) & (sum(corrPC_one) == sum(corrPC_two))) {
        
        finalHapMerged_Ordered <- list()
        for (hapNb in 1:length(indexOverlappingHap)) {
          
          #merging the full HQ with the lowHQ haplotype
          if (length(setdiff(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))) >= 1) {
            
            lastVarCommon <- intersect(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))[1]
            varOut <- intersect(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))[-1]
            index_varOut <- match(varOut, colnames(hapList[[indexOverlappingHap[hapNb]]]))
            
            finalHapMerged <- merge(newHap[[hapNb]], hapList[[indexOverlappingHap[hapNb]]][,-index_varOut], by = lastVarCommon)
            finalHapMerged_Ordered[[hapNb]] <- finalHapMerged[,order(as.numeric(colnames(finalHapMerged)), decreasing = F)]
            
          } else if (length(setdiff(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))) == 0) {
            
            finalHapMerged_Ordered[[hapNb]] <- newHap[[hapNb]]
          }
          #return(finalHapMerged_Ordered)
        }
        
        
        cat(paste("Variant:", LQVarCoord, "kept.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
        if (length(finalHapMerged_Ordered) == 1) {
          
          finalHap <- unique(finalHapMerged_Ordered)
          
        } else if (length(finalHapMerged_Ordered) > 1) {
          
          finalHap <- excludingDuplicates(unique(finalHapMerged_Ordered))[[1]]
          
        }
        
        return(finalHap)
        
      } else {
        
        cat(paste("Variant:", LQVarCoord, "not supported or ERROR: haplotypes differ.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
        finalHap <- 0
        return(finalHap)
      }
    } else if (LQhapQC == "restricted") {
      
      cat(paste("Variant:", LQVarCoord, "not supported or ERROR: less than five comparisons.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      finalHap <- 0
      return(finalHap)
    } else {
      
      cat(paste("Variant:", LQVarCoord, "not supported or ERROR: haplotypes differ.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      finalHap <- 0
      return(finalHap)
      
    }
    
  } else {
    
    cat(paste("Variant:", LQVarCoord, "not supported. ERROR: not enough HQ haplotypes", sep=" "), file=outputLQPhasing, sep="\n", append=T)
    finalHap <- 0
    return(finalHap)
    
    
  }
  
}



#######################
## NEW FUNCTIONS - Aug 18 2017
##
#######################
#this function takes the subMatrix following QC and generates the possible haplotypes between a provided LQ variant and all the possible
#and previously built HQ haplotypes. It exports a list with the haplotype
#it takes global variables: hapList, subMatrix, LQhapQC
creatingLQhapsWithTags <- function(LQVarCoord) {
  
  #pos <- 3
  #LQVarCoord <- varLQstep[2]
  #checks how many pairwise comb in the LQ matrix are related to HQ haplotypes and for how many sites within each HQ hap are covered by the LQ list
  #overlappingHap is a vector of the length of the hapList and the numbers represent the nb of sites with links information in the LQ matrix
  overlappingHap <- unlist(lapply(1:length(hapList), function(x) { sum(is.element(subMatrix[which(subMatrix[,2] == LQVarCoord),3], colnames(hapList[[x]]))) }))
  indexOverlappingHap <-  which(overlappingHap >= 2) #note before was >=
  #print(length(indexOverlappingHap))
  
  
  if (length(indexOverlappingHap) > 0) {
    
    corrPC_one <- c()
    corrPC_two <- c()
    newHap <- list()  
    #hapNb <- 75
    countHapOver <- 1
    for (hapNb in indexOverlappingHap) {
      
      #subsetting the subMatrix per LQ variant
      perHQhap_NewLQPos_m <- subMatrix[intersect(which(subMatrix[,2] == LQVarCoord), which(subMatrix[,3] %in% as.numeric(colnames(hapList[[hapNb]])))),]
      #subsetting the matrix with the links' order
      perHQhap_NewLQPos_order_m <- subOrder[intersect(which(subMatrix[,2] == LQVarCoord),which(subMatrix[,3] %in% as.numeric(colnames(hapList[[hapNb]])))),]
      vectorList <- list()
      
      
      for (col in 1:nrow(perHQhap_NewLQPos_m)) {
        
        linkSupported <- perHQhap_NewLQPos_order_m[col,][1:2]
        
        
        
        if (is.element(1, linkSupported) & is.element(4, linkSupported)) {
          
          vectorList[[col]] <- rbind(c(0,0),c(1,1))
          colnames(vectorList[[col]]) <- c(as.character(perHQhap_NewLQPos_m[col,2]), as.character(perHQhap_NewLQPos_m[col,3]))
          
        } else if (is.element(2, linkSupported) & is.element(3, linkSupported)) {
          
          vectorList[[col]] <- rbind(c(0,1),c(1,0))
          colnames(vectorList[[col]]) <- c(as.character(perHQhap_NewLQPos_m[col,2]), as.character(perHQhap_NewLQPos_m[col,3]))
          
        }
        
        if (col == 1) {
          
          newHap[[countHapOver]] <- vectorList[[col]]
          
        } else {
          
          newHap[[countHapOver]] <- merge(newHap[[countHapOver]], vectorList[[col]], by = LQVarCoord)
          
          
        }
        
      }
      
      # the first which does this - looking for the row number with a zero in the first overlapping site between the HQ hap and the new LQ hap that is not the LQ SITE
      #the following vector contain the combination of states followed by the first zero in the HQ haplotype
      hapZero_One <- as.vector(hapList[[hapNb]][which(hapList[[hapNb]][,colnames(hapList[[hapNb]])[colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])][1]] == 0 ),colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])])
      #the following vector contain the combination of states followed by the first zero in the LQ haplotype
      hapZero_Two <- as.vector(newHap[[countHapOver]][which(newHap[[countHapOver]][,colnames(newHap[[countHapOver]])[colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])][1]] == 0 ),colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])])
      #example:
      #HQ hap 
      # SNP1 SNP2
      # 0    1
      # 1    0
      #SNP1 and 2 are overlapping between the HQ hap and the new LQ hap
      # hapZero_One contains coordinate c(0, column index(SNP1, SNP2))
      
      #the two vectors below contain the same as before but the combination starting with 1.
      hapOne_One <- as.vector(hapList[[hapNb]][which(hapList[[hapNb]][,colnames(hapList[[hapNb]])[colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])][1]] == 1 ),colnames(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])])
      hapOne_Two <- as.vector(newHap[[countHapOver]][which(newHap[[countHapOver]][,colnames(newHap[[countHapOver]])[colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])][1]] == 1 ),colnames(newHap[[countHapOver]]) %in% colnames(hapList[[hapNb]])])
      
      #here we compare the HQ hap with the LQ hap for the overlapping sites
      corrPC_one[countHapOver] <- sum(apply(m <- rbind(hapZero_One, hapZero_Two), 2, function(x) { is.element(x[1], x[2])}))
      corrPC_two[countHapOver] <- sum(apply(m <- rbind(hapOne_One, hapOne_Two), 2, function(x) { is.element(x[1], x[2])}))
      countHapOver <- countHapOver+1
      
    } #close for across compatible HQ haps
    
    if (sum(corrPC_one) == sum(corrPC_two)) { #added Aug 17
      
      elementIDToKeep <- which(corrPC_one == overlappingHap[indexOverlappingHap])
      elementIDToRm <- which(corrPC_one != overlappingHap[indexOverlappingHap])
      if (length(elementIDToRm) > 0) {
        newHap <- newHap[-c(elementIDToRm)]
      }
      #corrPC_one/two contain the nb of EQUAL comparisons between the new and the HQ haplotype
      #if the nb of comparisons is >= 5 the LQ is kept and the LQ haploytpe built
      if (length(elementIDToKeep) > 0) { #added Aug 17
        
        indexOverlappingHap <- indexOverlappingHap[elementIDToKeep]
        #if ((sum(corrPC_one)[elementIDToKeep]/sum(overlappingHap[indexOverlappingHap]) == 1)) {
        if (sum(overlappingHap[indexOverlappingHap]) >= 5) {
          
          finalHapMerged_Ordered <- list()
          for (hapNb in 1:length(indexOverlappingHap)) {
            
            #merging the full HQ with the lowHQ haplotype
            if (length(setdiff(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))) >= 1) {
              
              lastVarCommon <- intersect(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))[1]
              varOut <- intersect(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))[-1]
              index_varOut <- match(varOut, colnames(hapList[[indexOverlappingHap[hapNb]]]))
              
              finalHapMerged <- merge(newHap[[hapNb]], hapList[[indexOverlappingHap[hapNb]]][,-index_varOut], by = lastVarCommon)
              finalHapMerged_Ordered[[hapNb]] <- finalHapMerged[,order(as.numeric(colnames(finalHapMerged)), decreasing = F)]
              
            } else if (length(setdiff(colnames(hapList[[indexOverlappingHap[hapNb]]]),colnames(newHap[[hapNb]]))) == 0) {
              
              finalHapMerged_Ordered[[hapNb]] <- newHap[[hapNb]]
            }
            #return(finalHapMerged_Ordered)
          }
          
          
          listHapTags <- lapply(1:length(indexOverlappingHap), function(x) { taggingLQvar_HQhap_combination(LQVarCoord, finalHapMerged_Ordered, x, indexOverlappingHap) })
          listDoubleHaps <- list()
          if (length(indexOverlappingHap) == 1) {
            
            tmpTag <- unlist(strsplit(listHapTags[[1]], split="_"))
            listDoubleHaps <- list(paste(tmpTag[1], tmpTag[2], "X", tmpTag[3], "X", sep = "_")) 
            vListDoubleHaps <- unlist(listDoubleHaps)
            return(vListDoubleHaps)
            
          } else if (length(indexOverlappingHap) >= 2) {
            
            df_hapTags <- data.frame(matrix(unlist(lapply(1:length(listHapTags), function(x) {strsplit(listHapTags[[x]], split = "_")})), ncol = 3, byrow = T))
            
            countRow <- 1
            
            while(countRow < nrow(df_hapTags)) {
              
              previousHapTmp <- df_hapTags[countRow,2]
              nextHapTmp <- df_hapTags[countRow+1,2]
              firstComb <- df_hapTags[countRow,3]
              
              if (df_hapTags[countRow,3] == "AD" & df_hapTags[countRow+1,3] == "AD") {
                
                secondComb <- "AD"
              } else if (df_hapTags[countRow,3] == "AD" & df_hapTags[countRow+1,3] == "BC") {
                secondComb <- "BC"
                
              } else if (df_hapTags[countRow,3] == "BC" & df_hapTags[countRow+1,3] == "AD") {
                secondComb <- "BC"
                
              } else if (df_hapTags[countRow,3] == "BC" & df_hapTags[countRow+1,3] == "BC") {
                secondComb <- "AD"
              }
              listDoubleHaps[[countRow]] <- paste(df_hapTags[1,1], previousHapTmp, nextHapTmp, firstComb, secondComb, sep="_")
              countRow <- countRow+1
            }
            
            cat(paste("Variant:", LQVarCoord, "kept.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
            
            vListDoubleHaps <- unlist(listDoubleHaps)
            return(vListDoubleHaps)
          }
        } else {
          
          cat(paste("Variant:", LQVarCoord, "not supported or ERROR: less than five compatible comparisons.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
          vListDoubleHaps <- 0
          return(vListDoubleHaps)
        }
        
      } else {
        
        cat(paste("Variant:", LQVarCoord, "not supported or ERROR: no HQ hap with compatible comparisons.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
        vListDoubleHaps <- 0
        return(vListDoubleHaps)
        
      }
    } else {
      
      cat(paste("Variant:", LQVarCoord, "not supported. ERROR: HQ haps are wrong.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      vListDoubleHaps <- 0
      return(vListDoubleHaps)
      
    }
    
  } else {
    
    cat(paste("Variant:", LQVarCoord, "not supported. ERROR: not enough HQ haplotypes", sep=" "), file=outputLQPhasing, sep="\n", append=T)
    vListDoubleHaps <- 0
    return(vListDoubleHaps)
    
  }
  
}
#####--------------
##---------------

taggingLQvar_HQhap_combination <- function(LQvar, list_full_HQhaps, nb_of_haplotype, indexHQhap) { #takes global variable --> hapList 
  
  # LQvar <- LQVarCoord
  # list_full_HQhaps <- finalHapMerged_Ordered
  # nb_of_haplotype <-  1 
  # indexHQhap <-  indexOverlappingHap
  
  LQvarIndex <- which(names(list_full_HQhaps[[nb_of_haplotype]]) == LQvar)
  rowIndex_Zero_in_LQvarIndex <- which(list_full_HQhaps[[nb_of_haplotype]][,LQvarIndex] == 0)
  
  firstSiteHQhapindex <- which(names(list_full_HQhaps[[nb_of_haplotype]]) == colnames(hapList[[indexHQhap[[nb_of_haplotype]]]])[1])
  if (length(firstSiteHQhapindex) == 0) {
    print("Error: HapList names not found.")
    break;
  }
  
  hapTmp <- list_full_HQhaps[[nb_of_haplotype]][rowIndex_Zero_in_LQvarIndex, c(LQvarIndex,firstSiteHQhapindex)]
  hapID <- c() 
  if (sum(hapTmp) == 0) { #the LQvar is always the REF, meaning is assumed to be zero
    
    hapID <- "AD"
  } else if (sum(hapTmp) == 1) {
    hapID <- "BC"
  }
  
  hapCombination <- paste(LQvar, indexHQhap[[nb_of_haplotype]], hapID, sep="_")
  return(hapCombination)
}

#######################
##
##
#######################

buildingHQ_connections <- function(hqListSize) { #takes hapList and hap_df as variables.
  
  
  HQhap_index_x <- 1
  HQhap_index_y <- HQhap_index_x+1
  linksTags <- c("AD", "BC", "X")
  tmpConnTag <- c()
  
  while (HQhap_index_x < hqListSize) {
    
    tmp_v <- as.vector(table(hap_df[which(hap_df[,2] == HQhap_index_x & hap_df[,3] == HQhap_index_y),5]))
    
    if (sum(tmp_v) > 5 & tmp_v[which(tmp_v == max(tmp_v))] >= 5*sum(tmp_v[-which(tmp_v == max(tmp_v))])) {
      
      
      tmpConnTag <- c(tmpConnTag, paste(HQhap_index_x,  HQhap_index_y, linksTags[which(tmp_v == max(tmp_v))], sep = "_"))
      HQhap_index_x <- HQhap_index_y
      HQhap_index_y <- HQhap_index_x+1
      
    } else {
      print("Error: not enough comparisons.")
      HQhap_index_y <- HQhap_index_y+1
    }
    if (HQhap_index_y > hqListSize) {
      break;
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
  
  #hapNb <- 1
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
  
}

#################################
##
## END OF FUNCTIONS
##
#################################
