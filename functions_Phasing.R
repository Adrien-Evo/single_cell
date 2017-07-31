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

#goes around the hap list until they are all merged and not more sites belong to two alternative chunks
excludingDuplicates <- function(tmpUniqueHapList) {
  
  #tmpUniqueHapList <- unique_h_list
  #tmpUniqueHapList <- newHap
  countDupList <- 1
  
  
while(T) { # haplopytes along the list
    
    tmp_hapList <- mergingHaplotypes(tmpUniqueHapList)
    tmp_namesHapList <- namesHapListFunction(tmp_hapList)
    
    if (sum(duplicated(tmp_namesHapList)) == 0) {
      
      print("Merging step, done!")
      break;
      
    } else {
      
      print(paste("Overlaping list nb:", countDupList, sep=" "))
      tmpUniqueHapList <- tmp_hapList
      countDupList <- countDupList+1
    }
    
    
  }
  return(list(tmp_hapList, tmp_namesHapList))
}



##-----------
#--------
subsettingByQuality_computingHQLinks <- function(genotypeMatrix) { #it takes the geno M filtered for those variants with > minHQCells 
  #takes the GLOBAL variables: minNbCells, output, linksOutput
  #it does not return anything, it rather creates a file with all links following the qual. conditions
  
  HQ_genotypeMatrix <- genotypeMatrix
  count_var_x <- 1
  count_var_y <- 2
  nbWindows <- 0
  
  
  while (count_var_x < nrow(HQ_genotypeMatrix)) {
    
    count_var_y <- count_var_x+1
    while (count_var_y <= nrow(HQ_genotypeMatrix)) {
      
      tmp <- HQ_genotypeMatrix[c(count_var_x,count_var_y),]
      col_i <- c()
      indexToTake <- as.vector(which(apply(tmp, 2, function(x) { col_i <- c(col_i, sum(is.na(x))) }) == 0))
      print(paste("Nb of overlapping cells:", length(indexToTake)))
      if (length(indexToTake) >= minNbCells) {
        
        tmp <- tmp[,indexToTake]
        v_countLinks <- apply(vapply(1:ncol(tmp), FUN = get_LinkCount, FUN.VALUE = rep(0,4)), 1, sum)
        print(v_countLinks)
        
        l <- order(v_countLinks, decreasing = T)
        
        ratio <- v_countLinks[l[2]]/5
        
        #if (v_countLinks[l[1]] >= minNbLinks & v_countLinks[l[2]] >= minNbLinks & v_countLinks[l[3]] < ratio) {
        
        cat("Link found", file=output, sep="\n", append = T)
        nbWindows <- nbWindows+1
        cat(paste("Link nb:", nbWindows, sep=" "), file=output, sep="\n", append=T)
        cat(paste("Two SNP window:", varNames[count_var_x], "|", varNames[count_var_y], sep=" "), file=output, sep="\n", append = T)
        
        cat(paste("chr15", varNames[count_var_x], varNames[count_var_y], paste("0/0/", v_countLinks[1], sep =""), paste("0/1/", v_countLinks[2], sep =""), paste("1/0/", v_countLinks[3], sep =""), paste("1/1/", v_countLinks[4], sep =""), collapse ="\t"), file=linksOutput, sep="\n", append=T)
      }
      count_var_y <- count_var_y+1
    }
    count_var_x <- count_var_x+1
    
  }
  
  
}
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
##--------------------
#--------------

#creating a matrix with all overlapping haps 
creatingHapMatrixFromLinkCountMatrix <- function(linkCountM_subset, orderLinkCountM_subset) { #it takes subMatrix and subOrder and return a matrix with all haplotypes
  
  nrow(linkCountM_subset)
  names_subM <- sort(union(linkCountM_subset[,1],linkCountM_subset[,2]))
  
  for (col in 1:nrow(linkCountM_subset)) {
    
    focalPosFirst <- linkCountM_subset[col,1]
    focalPosSec <- linkCountM_subset[col,2]
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

subsettingByQuality_computingLQLinks <- function(outFName) {
  
  
  cat(paste("PosOne", "PosTwo", "0/0", "0/1", "1/0", "1/1", sep="\t"), file=outFName, sep="")
  varNb <- 1
  
  while(varNb<length(hetSNPsPos)) { #######CHANGE HERE
    
    #varNb <- 1
    
    if (!is.element(hetSNPsPos[varNb], colnames(matrixSortedHap))) {
      
      
      subGenoQ <- openGenoq[match(c(hetSNPsPos[varNb],colnames(matrixSortedHap)), varNames),]
      subGeno <- openGeno[match(c(hetSNPsPos[varNb],colnames(matrixSortedHap)), varNames),]
      indexHQ <- lapply(1:nrow(subGenoQ), function(x) { which(subGenoQ[x,] > 20 )})
      
      df_countLinks <- list()
      countComparisons <- 2
      countLinks <- 0
      
      while(countComparisons <= nrow(subGeno)) {
        
        print(countComparisons)
        interSectVar <- intersect(indexHQ[[1]], indexHQ[[countComparisons]])
        if (length(interSectVar) >= 4) {
          
          countLinks <- countLinks+1     
          tmp <- subGeno[c(1,countComparisons),interSectVar]
          tmpLinks <- apply(vapply(1:ncol(tmp), FUN = get_LinkCount, FUN.VALUE = rep(0,4)), 1, sum)
          
          df_countLinks[[countLinks]] <- matrix(cbind(as.numeric(hetSNPsPos[varNb]),as.numeric(colnames(matrixSortedHap)[(countComparisons-1)]),matrix(tmpLinks, ncol=4)), ncol=6)
        }
        countComparisons <- countComparisons+1
        
      }
      
      write.table(matrix(unlist(df_countLinks), ncol=6, byrow = T), file=outFName, append = T, quote = F, row.names = F, col.names = F, sep="\t")
      varNb <- varNb+1
    } else {
      
      varNb <- varNb+1
      
      
    } 
  }
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
  indexOverlappingHap <-  which(overlappingHap >= 2) #note before was >=
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

#################################
##
## END OF FUNCTIONS
##
#################################
