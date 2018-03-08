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
  overlappingHap <- unlist(lapply(hapList, function(hapl) { sum(is.element(list_counts_order[[1]][,3], names(hapl))) })) # modified by Nov 27
  indexOverlappingHap <-  which(overlappingHap >= 2) #note before was >=
  #print(length(indexOverlappingHap))
  #getting the HQ haps closer from the LQ variant
  #which(order(sapply(1:length(hapList), function(x) { abs(as.numeric(LQVarCoord)-as.numeric(colnames(hapList[[x]])[1])) } )) <= 10)
  closerHaps <- order(sapply(1:length(hapList), function(x) { min(abs(as.numeric(LQVarCoord)-as.numeric(names(hapList[[x]])))) } ))[1:nbOfCloserHQhap] #modified by Jan 8
  intersectCov_closerHaps <- intersect(closerHaps, indexOverlappingHap)
  LQ_HQ_tags <- rep(NA, length(overlappingHap))
  
  overlappingHap[intersectCov_closerHaps][which(overlappingHap[intersectCov_closerHaps] > 5)] <- 5
  #print(intersectCov_closerHaps)
  
  if (length(intersectCov_closerHaps) > 0 & sum(overlappingHap[intersectCov_closerHaps]) >= 5) {
    
    corrPC_one <- c()
    corrPC_two <- c()
    newHap <- list()  
    #hapNb <- 75
    countHapOver <- 1
    for (hapNb in intersectCov_closerHaps) {
      
      #subsetting the subMatrix per LQ variant
      perHQhap_NewLQPos_m <- list_counts_order[[1]][which(list_counts_order[[1]][,3] %in% as.numeric(names(hapList[[hapNb]]))),] #modified by late Nov 27
      #subsetting the matrix with the links' order
      perHQhap_NewLQPos_order_m <- list_counts_order[[2]][which(list_counts_order[[1]][,3] %in% as.numeric(names(hapList[[hapNb]]))),] #modified by late Nov 27
      
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
      
      newHap[[countHapOver]] <-  newHap[[countHapOver]][c(which(newHap[[countHapOver]][,1] == 0), which(newHap[[countHapOver]][,1] == 1)),]
      #TO TEST TOMORROW - Jan 11
      #overlappingHap[intersectCov_closerHaps][which(overlappingHap[intersectCov_closerHaps] > 5)] <- 5
      
      
      
      # the first which does this - looking for the row number with a zero in the first overlapping site between the HQ hap and the new LQ hap that is not the LQ SITE
      #the following vector contain the combination of states followed by the first zero in the HQ haplotype
      HQHap_zero <- hapList[[hapNb]][names(hapList[[hapNb]]) %in% colnames(newHap[[countHapOver]])]
      #the following vector contain the combination of states followed by the first zero in the LQ haplotype
      LQ_HQhap_zero <- as.vector(newHap[[countHapOver]][1, colnames(newHap[[countHapOver]]) %in% names(hapList[[hapNb]])])
      #example:
      
      similarity <- sum(apply(m <- rbind(HQHap_zero, LQ_HQhap_zero), 2, function(x) { is.element(x[1], x[2])}))/overlappingHap[hapNb]
      
      if (similarity < 0.20) {
        
        tag <- -1
        
      } else if (similarity >= 0.80) {
        
        tag <- 1
        
      } else {
        
        tag <- NA
        
      }
      
      LQ_HQ_tags[hapNb] <- tag
      countHapOver <- countHapOver+1
      
    } #close for across compatible HQ haps
    
    
    if (sum(is.na(LQ_HQ_tags))/length(LQ_HQ_tags) < 1) { #CHANGE TOMORROW TO >= 0.8 added by Jan17
      
      
      cat(paste("Variant:", LQVarCoord, "kept.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      #generating random file
      tmpfile <- paste(sample(c(LETTERS,letters, 0:9), 20, replace=TRUE), collapse="") #changed by March 2
      cat(paste(LQVarCoord, paste0(LQ_HQ_tags, collapse = "\t"), sep = "\t"), file=paste0(paste(pathToTmp, tmpfile, sep = "/"), ".txt"), sep="\n", append = T) #changed by March 2


    } else {
      
      cat(paste("Variant:", LQVarCoord, "not supported or ERROR: hap length < 2 and total sites < 5.", sep=" "), file=outputLQPhasing, sep="\n", append=T)
      # finalHap <- 0
      # return(finalHap)
    }
    
  } else {
    
    cat(paste("Variant:", LQVarCoord, "not supported. ERROR: not enough HQ haplotypes", sep=" "), file=outputLQPhasing, sep="\n", append=T)
    # finalHap <- 0
    # return(finalHap)
    
  }
  
}