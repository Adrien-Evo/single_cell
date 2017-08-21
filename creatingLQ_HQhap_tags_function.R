#######################
## NEW FUNCTIONS - Aug 18 2017
##
#######################
#this function takes the subMatrix following QC and generates the possible haplotypes between a provided LQ variant and all the possible
#and previously built HQ haplotypes. It exports a list with the haplotype
#it takes global variables: hapList, subMatrix, subOrder
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