#######################
## NEW FUNCTIONS - july 27 2017
##
#######################
#this function takes the subMatrix following QC and generates the possible haplotypes between a provided LQ variant and all the possible
#and previously built HQ haplotypes. It exports a list with the haplotype
#it takes global variables: hapList, subMatrix, nbOfCloserHQhap
creatingLQhaplotypes <- function(list_counts_order) {
  
  #pos <- 3
  #list_counts_order <- links_and_order_List[c(545)]
  LQVarCoord <- as.character(unique(list_counts_order[[1]][,2])) # modified by Nov 27
  cat(paste("Analysing variant:", LQVarCoord, ".", sep=" "), file=outputLQPhasing, sep="\n", append=T)
  #LQVarCoord <- "25032267"
  #checks how many pairwise comb in the LQ matrix are related to HQ haplotypes and for how many sites within each HQ hap are covered by the LQ list
  #overlappingHap is a vector of the length of the hapList and the numbers represent the nb of sites with links information in the LQ matrix
  overlappingHap <- unlist(lapply(hapList, function(hapl) { sum(is.element(list_counts_order[[1]][,3], colnames(hapl))) })) # modified by Nov 27
  indexOverlappingHap <-  which(overlappingHap >= 2) #note before was >=
  #print(length(indexOverlappingHap))
  #getting the HQ haps closer from the LQ variant
  #which(order(sapply(1:length(hapList), function(x) { abs(as.numeric(LQVarCoord)-as.numeric(colnames(hapList[[x]])[1])) } )) <= 10)
  closerHaps <- order(sapply(1:length(hapList), function(x) { min(abs(as.numeric(LQVarCoord)-as.numeric(colnames(hapList[[x]])))) } ))[1:nbOfCloserHQhap] #modified by Jan 8
  intersectCov_closerHaps <- intersect(closerHaps, indexOverlappingHap)
  #print(intersectCov_closerHaps)
  
  if (length(intersectCov_closerHaps) > 0 & sum(overlappingHap[intersectCov_closerHaps]) >= 5) {
    
    corrPC_one <- c()
    corrPC_two <- c()
    newHap <- list()  
    #hapNb <- 75
    countHapOver <- 1
    for (hapNb in intersectCov_closerHaps) {
      
      #subsetting the subMatrix per LQ variant
      perHQhap_NewLQPos_m <- list_counts_order[[1]][which(list_counts_order[[1]][,3] %in% as.numeric(colnames(hapList[[hapNb]]))),] #modified by late Nov 27
      #subsetting the matrix with the links' order
      perHQhap_NewLQPos_order_m <- list_counts_order[[2]][which(list_counts_order[[1]][,3] %in% as.numeric(colnames(hapList[[hapNb]]))),] #modified by late Nov 27
      
      #TO TEST TOMORROW  - Jan 11
      if (nrow(perHQhap_NewLQPos_m) > 5) {
        
        perHQhap_NewLQPos_m_tmp <- perHQhap_NewLQPos_m[order(abs(as.numeric(LQVarCoord)-as.numeric(as.character(perHQhap_NewLQPos_m[,3]))), decreasing = F)[1:5],] #added by Jan 8
        perHQhap_NewLQPos_order_m_tmp <- perHQhap_NewLQPos_order_m[order(abs(as.numeric(LQVarCoord)-as.numeric(as.character(perHQhap_NewLQPos_m[,3]))), decreasing = F)[1:5],]  #added by Jan 8
        perHQhap_NewLQPos_m <- perHQhap_NewLQPos_m_tmp[order(as.numeric(as.character(perHQhap_NewLQPos_m_tmp[,3]))),] #added by Jan 9
        perHQhap_NewLQPos_order_m <- perHQhap_NewLQPos_order_m_tmp[order(as.numeric(as.character(perHQhap_NewLQPos_m_tmp[,3]))),] #added by Jan 9
      }
      
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
          
          newHap[[countHapOver]] <- merge(newHap[[countHapOver]], vectorList[[col]], by = as.character(LQVarCoord))
          
        }
        
      }
      
      #TO TEST TOMORROW - Jan 11
      overlappingHap[intersectCov_closerHaps][which(overlappingHap[intersectCov_closerHaps] > 5)] <- 5
      
      
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
      
      compatibleComp_Zero_Zero <- sum(apply(m <- rbind(hapZero_One, hapZero_Two), 2, function(x) { is.element(x[1], x[2])}))
      compatibleComp_Zero_One <- sum(apply(m <- rbind(hapZero_One, hapOne_Two), 2, function(x) { is.element(x[1], x[2])}))
      
      #here we compare the HQ hap with the LQ hap for the overlapping sites
      corrPC_one[countHapOver] <-  max(compatibleComp_Zero_Zero, compatibleComp_Zero_One)
      
      
      if (which(c(compatibleComp_Zero_Zero, compatibleComp_Zero_One) == corrPC_one[countHapOver]) == 1) {
        
        indexToKeepNewHap <- as.vector(which(apply(m <- rbind(hapZero_One, hapZero_Two), 2, function(x) { is.element(x[1], x[2])}) == TRUE))+1
        
      } else if (which(c(compatibleComp_Zero_Zero, compatibleComp_Zero_One) == corrPC_one[countHapOver]) == 2) {
        
        indexToKeepNewHap <- as.vector(which(apply(m <- rbind(hapZero_One, hapOne_Two), 2, function(x) { is.element(x[1], x[2])}) == TRUE))+1
      }
    
      newHap[[countHapOver]] <- newHap[[countHapOver]][,c(1,indexToKeepNewHap)] #added by Jan 15
      #corrPC_two[countHapOver] <- sum(apply(m <- rbind(hapOne_One, hapOne_Two), 2, function(x) { is.element(x[1], x[2])}))
      
      countHapOver <- countHapOver+1
      
    } #close for across compatible HQ haps
    
    #corrPC_one/two contain the nb of EQUAL comparisons between the new and the HQ haplotype
    #if the nb of comparisons is >= 5 the LQ is kept and the LQ haploytpe built
    #if () { #modified by late Nov 27
    newHap <- newHap[which(corrPC_one/overlappingHap[intersectCov_closerHaps] >= 0.80)] #added Jan 15
    intersectCov_closerHaps <- intersectCov_closerHaps[which(corrPC_one/overlappingHap[intersectCov_closerHaps] >= 0.80)] #added Jan 17
    
    if (length(newHap) > 0) { #CHANGE TOMORROW TO >= 0.8 added by Jan17
      
      finalHapMerged_Ordered <- list()
      for (hapNb in 1:length(intersectCov_closerHaps)) {
        
        #merging the full HQ with the lowHQ haplotype
        if (length(setdiff(colnames(hapList[[intersectCov_closerHaps[hapNb]]]),colnames(newHap[[hapNb]]))) >= 1) {
          
          lastVarCommon <- intersect(colnames(hapList[[intersectCov_closerHaps[hapNb]]]),colnames(newHap[[hapNb]]))[1]
          varOut <- intersect(colnames(hapList[[intersectCov_closerHaps[hapNb]]]),colnames(newHap[[hapNb]]))[-1]
          index_varOut <- match(varOut, colnames(hapList[[intersectCov_closerHaps[hapNb]]]))
          
          
          #cat(paste("Length of index to rm:", length(index_varOut), sep=" "), file=outputLQPhasing, sep="\n", append=T) #add Jan 15
          finalHapMerged <- merge(newHap[[hapNb]], hapList[[intersectCov_closerHaps[hapNb]]][,-index_varOut], by = lastVarCommon)
          finalHapMerged_Ordered[[hapNb]] <- finalHapMerged[,order(as.numeric(colnames(finalHapMerged)), decreasing = F)]
          
        } else if (length(setdiff(colnames(hapList[[intersectCov_closerHaps[hapNb]]]),colnames(newHap[[hapNb]]))) == 0) {
          
          finalHapMerged_Ordered[[hapNb]] <- newHap[[hapNb]]
        }
        #return(finalHapMerged_Ordered)
      }
      
      cat(paste("Variant:", LQVarCoord, "kept.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      if (length(finalHapMerged_Ordered) == 1) {
        
        finalHap <- unique(finalHapMerged_Ordered)[[1]] #changed [[1]] Oct 4
        
      } else if (length(finalHapMerged_Ordered) > 1) {
        
        finalHap <- excludingDuplicates(unique(finalHapMerged_Ordered))[[1]]
        
      }
      
      return(finalHap)
      
    } else {
      
      cat(paste("Variant:", LQVarCoord, "not supported or ERROR: hap length < 2 and total sites < 5.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      # finalHap <- 0
      # return(finalHap)
    }
    # } else {
    #   
    #   cat(paste("Variant:", LQVarCoord, "not supported or ERROR: less than five comparisons.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
    #   # finalHap <- 0
    #   # return(finalHap)
    # } 
    
  } else {
    
    cat(paste("Variant:", LQVarCoord, "not supported. ERROR: not enough HQ haplotypes", sep=" "), file=outputLQPhasing, sep="\n", append=T)
    # finalHap <- 0
    # return(finalHap)
    
  }
  
}
