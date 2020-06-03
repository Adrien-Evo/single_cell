##################################
##################################
##
##          Functions
##
##################################
##################################

remove_errors_from_list_match <- function(i) { #global: list_match, seqList
  #i <- 93
  match_v <- list_match[[i]]
  if(length(match_v) > 1) {
    if(sum(match_v == TRUE) > sum(match_v == FALSE)) {
      status <- TRUE
      
    } else {
      status <- FALSE
    }
    
    index_dissimilar <- which(match_v != status)
    if(length(seqList[[(i-1)]]) == length(index_dissimilar)) { #if there are no 0 or 100% missmatches
      #or the missmatches detected are sequences and are not real errors.
      newMatch <- list_match[[i]]
      nb_errors <- 0
      
    } else if (length(seqList[[(i-1)]]) == 0 & length(seqList[[(i-1)]]) < length(index_dissimilar)) { 
      #The missmatches detected are ALL errors
      
      newMatch <- list_match[[i]][-index_dissimilar]
      nb_errors <- length(index_dissimilar)
      
    } else if (length(seqList[[(i-1)]]) > 0 & length(seqList[[(i-1)]]) < length(index_dissimilar)) { #if there are sequences identified.
      
      index_errors <- setdiff(index_dissimilar, seqList[[(i-1)]])
      newMatch <- list_match[[i]][-index_errors]
      nb_errors <- length(index_errors)
      
    }
    return(list(cleanMatches=newMatch, errorNb=nb_errors))
    
  } else if(length(match_v) == 1 & is.na(match_v)) {
    return(list(cleanMatches=NA, errorNb=NA))
  }
}


#function to identify sequences in the vectors 
#comparing one target cell with all the others
#this function removed sequeces of length = 2 
#where the two SNPs fall within a distance of 150bp
# (likely in the same read)
retrieve_sequences <- function(match_v) { #global: NONE
  
  if(length(match_v) > 1) { #otherwise it is NA and there was no comparison
    #match_v <- list_match[[22]]
    if(sum(match_v == TRUE) > sum(match_v == FALSE)) {
      status <- TRUE
      
    } else {
      status <- FALSE
    }
    
    index_dissimilar <- which(match_v != status)
    start <- 1
    nextPos <- start+1
    seqIdx <- list()
    COUNTSEQ <- 0
    seqStatus <- "no"
    
    while(start < length(index_dissimilar)){
      
      #print(paste0("Starting with: ", index_dissimilar[start]))
      
      if(index_dissimilar[nextPos] == (index_dissimilar[start]+1)) {
        if (seqStatus == "no") {
          COUNTSEQ <- COUNTSEQ+1
          seqIdx[[COUNTSEQ]] <- c(index_dissimilar[start], index_dissimilar[nextPos])
        } else {
          seqIdx[[COUNTSEQ]] <- c(seqIdx[[COUNTSEQ]], index_dissimilar[start], index_dissimilar[nextPos])
          seqIdx[[COUNTSEQ]] <- seqIdx[[COUNTSEQ]]
        }
        seqStatus <- "yes"
        
      } else {
        seqStatus <- "no"
        
      }
      start <- nextPos
      nextPos <- start+1    
    }
    
    #remove the seqs with two var at distances less then 150bp
    tmp_list <- lapply(seqIdx, function(seq_x) {
      
      #seq_x <- seqIdx[[1]]
      x <- unique(seq_x)
      names(x) <- unique(names(seq_x))
      
      if(length(x) == 2) {
        coord <- as.numeric(as.character(names(x)))
        gapSize <- coord[2]-coord[1]
        if(gapSize <= 150) {
          return(NULL)
        } else {
          return(x)
        }
        
      } else {
        return(x)
      }
    })
  #return(unlist(lapply(seqIdx, FUN = unique)))
  return(list(nb_of_seq=length(lengths(tmp_list)[which(lengths(tmp_list) != 0)]),seq=unlist(lapply(tmp_list, FUN = unique))))
  } else {
    return(list(nb_of_seq=NA, seq=NA))
  }
} 

#hmm function to infer parental states given the observation
inferring_bp_hmm <- function(cellIndex) {
  
  targetCell <- cellIndex
  #targetCell <- 135
  #subsetting the data
  indGeno <- droplevels(openGT[which(openGT[,targetCell+2] != "."), targetCell+2])
  positions <- openGT[which(openGT[,targetCell+2] != "."), 2]
  cellID <- colnames(openGT)[targetCell+2]
  print(paste0("Inferring recombination bp for cell: ", cellID))
  print("")
  inters_cell_parentHaps <- intersect(colnames(openParentalHaps), positions)
  subSetParentHaps <- openParentalHaps[,match(inters_cell_parentHaps, colnames(openParentalHaps))]
  subSetIndGeno <- indGeno[match(inters_cell_parentHaps, positions)]  
  subSetPos <- positions[match(inters_cell_parentHaps, positions)]
  
  #creating arrays
  #fullSeq is an array of the size of the phased haplotypes where the state path for each cell will be copied
  fullSeq <- rep(".", ncol(openParentalHaps))
  names(fullSeq) <- colnames(openParentalHaps)
  fullStateSeq <- rep(".", ncol(openParentalHaps))
  names(fullStateSeq) <- colnames(openParentalHaps)
  
  m_Seq <- matrix(rep(NA, (length(subSetIndGeno)+1)*2), nrow=2, ncol=length(subSetIndGeno)+1)
  colnames(m_Seq) <- c("start", as.character(subSetIndGeno))
  m_Seq[,1] <- c(0,0)
  path <- c()
  statePath <- c()
  
  for(l in 1:length(subSetIndGeno)) { #length(subSetIndGeno)
    #l <- 1
    error_tmp <- errorRates[which(errorRates[,1] == subSetPos[l]),2]
    if(is.na(error_tmp)) {
      error_tmp <- mean(errorRates[,2], na.rm = T)
    } 
    
    #i am considering that P is the upper haplotype
    stateP_upperHap <- subSetParentHaps[1,which(names(subSetParentHaps) == subSetPos[l])]
    stateM_lowerHap <- subSetParentHaps[2,which(names(subSetParentHaps) == subSetPos[l])]
    
    if(stateP_upperHap == "0") {
      #print("Paternal state: 0")
      emissionProb <- matrix(c(1-error_tmp, error_tmp, error_tmp, 1-error_tmp),nrow = 2, byrow = T)
    } else {
      #print("Paternal state: 1")
      emissionProb <- matrix(c(error_tmp, 1-error_tmp, 1-error_tmp, error_tmp),nrow = 2, byrow = T)
      
    }
    #emissionProb
    if(l > 1) {
      recombRat_tmp <- recombRate_perSite*(subSetPos[l]-subSetPos[l-1])
      state_trans_prob <- matrix(c(1-recombRat_tmp, recombRat_tmp, recombRat_tmp,1-recombRat_tmp),2)
    }
    
    if(l==1) {
      pM <- log2(emissionProb[which(states == "M"), which(subSetIndGeno[l]==symbols)])+log2(startProb)
      pP <- log2(emissionProb[which(states == "P"), which(subSetIndGeno[l]==symbols)])+log2(startProb)
    } else {
      pM <- log2(emissionProb[which(states == "M"), which(subSetIndGeno[l]==symbols)])+max(m_Seq[2,l]+log2(state_trans_prob[2,2]), m_Seq[1,l]+log2(state_trans_prob[1,2]))
      pP <- log2(emissionProb[which(states == "P"), which(subSetIndGeno[l]==symbols)])+max(m_Seq[1,l]+log2(state_trans_prob[1,1]), m_Seq[2,l]+log2(state_trans_prob[2,1]))
      
    }
    m_Seq[,l+1] <- c(pP, pM)
    if(length(which(c(pP,pM) == max(pP,pM))) == 2) {
      # print(l)
      # print(which(c(pP,pM) == max(pP,pM)))
      path <- c(path, sample(c(0,1), 1))
      statePath <- c(statePath, "U")
      #break;
    } else {
      path <- c(path,subSetParentHaps[which(c(pP,pM) == max(pP,pM)),l])      
      statePath <- c(statePath, states[which(c(pP,pM) == max(pP,pM))])
    }
    
  }
  
  names(path) <- subSetPos
  fullSeq[match(names(path), names(fullSeq))] <- path
  names(statePath) <- subSetPos
  fullStateSeq[match(names(statePath), names(fullStateSeq))] <- statePath
  
  SeqList <- list(cellName=cellID, inferredSeq=fullSeq, inferredStates=fullStateSeq)
  return(SeqList)
}

## function related to plotting inferred chunks
#################
#Function
#################
colPolygon <- function(coordTmp) {
  
  cTmp <- c()
  if (is.element(coordTmp, mCoord)) {
    
    cTmp <- "red"
  } else if (is.element(coordTmp, pCoord)) {
    cTmp <- "black"
  } else if (is.element(coordTmp, uCoord)) {
    cTmp <- "white"
  }
  return(cTmp)
}
#####---------

##Functions related to the computation of bp and parental chromosome chunks
compute_nb_chromosome_chunks <- function(cellIndx) {
  #Test function with one cell in vcf file
  #v <- cell_hmm_StateSeq_m[[1]
  #cellIndx <- 10
  print(cellIndx)
  v <- cell_hmm_StateSeq_m[[cellIndx]]
  #Stores first and last coordinate of each M/P block
  startM <- c()
  endM <- c()
  startP <- c()
  endP <- c()
  previousCoord <- c()
  #Neutral parental state of the chromosome (vs. M/P)
  state <- "I"
  COUNT_CHUNKS <- 0
  stateOrder <- c()
  
  #Scans cell from the first to last position of chromosome to detect parental origin (M or P)
  for (i in 1:length(v)) { #:length(v)
    #print(v[i])
    #Scans for first coordinate of M block
    if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores first coordinate of M block
      startM <- c(startM, as.numeric(names(v[i])))
      #print(paste0("Start of chrom chunk M: ", startM))
      #Scans for previous P block end coordinate up until new M block
      endP <- c(endP, previousCoord)
      stateOrder <- c(stateOrder,state)
      COUNT_CHUNKS <- COUNT_CHUNKS+1
      #Scans for first coordinate of new P block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores the first coordinate of P block 
      startP <- c(startP, as.numeric(names(v[i])))
      #print(paste0("Start of chrom chunk P: ", startP))
      COUNT_CHUNKS <- COUNT_CHUNKS+1
      stateOrder <- c(stateOrder,state)
      #Scans for previous M block end coordinate up until new P block
      endM <- c(endM, previousCoord)
    } 
    if (v[i] == "M" | v[i] == "P") {
      previousCoord <- as.numeric(names(v[i]))
      #print(previousCoord)
    }
    
  } #end of going along every element of the vector with the parental states
  
  lastCoord <-as.numeric(names(v[as.numeric(which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))])]))
  allCoord <- sort(c(startM, startP, endM, endP, lastCoord), decreasing = F)
  
  if(length(stateOrder) == (length(allCoord)/2)) {
    
    startCoord <- allCoord[seq(1,(length(allCoord)-1), by = 2)]
    endCoord <- allCoord[seq(2,length(allCoord), by = 2)]
    sizeChunks <- endCoord - startCoord 
    #print(paste0("Chunks detected. There are ", length(sizeChunks), " of sizes: ", sizeChunks))
  } else {
    print("Something went wrong...The number of coordinates is not even.")
  }
  return(list(chunksIdentity=stateOrder,sCoordinate=startCoord, eCoordinate=endCoord, chunkLength=sizeChunks))
}
####-----------------------
#this function removes parental chunks from the hmm list with fragment size smaller than the 
# fragLengthToExclude defined as global variable. 
remove_short_chunks <- function(chunkL_v, hmm_v, sCoor_v, eCoor_v) { #it takes fragLengthToExclude as global variable
  
  # ll <- 78
  # chunkL_v <- chunksLengths_clean[[ll]]
  # hmm_v <- list_hmm_clean[[ll]]
  # sCoor_v <- coordinatesStart_clean[[ll]]
  # eCoor_v <- coordinatesEnd_clean[[ll]]
  
  indx <- which(chunkL_v < fragLengthToExclude)
  if(length(indx) > 0) {
    if(is.element(1, indx)) {
      indx <- indx[-1]
    } else if(is.element(length(chunkL_v), indx)) {
      indx <- indx[-length(chunkL_v)]
    }
    if(length(indx) > 0){
      for(i in indx) {
        #stateSequences_clean[[which(chunksLengths_clean[[4]] < fragLengthToExclude)]]
        x0 <- sCoor_v[i]
        x1 <- eCoor_v[i]
        
        hmm_v[which(as.numeric(names(hmm_v)) >= x0 & as.numeric(names(hmm_v)) <= x1)] <- "."
        print(length(which(as.numeric(names(hmm_v)) >= x0 & as.numeric(names(hmm_v)) <= x1)))
      }
    }
  } 
  return(hmm_v)
}
####--------------
##count the number of SNPs inside non recombining chunks
counting_SNPs_per_chunk <- function(chunkL_v, hmm_v, sCoor_v, eCoor_v, cellID) { #it takes fragLengthToExclude as global variable
  
  # ll <- 78
  # chunkL_v <- chunksLengths_clean[[ll]]
  # hmm_v <- list_hmm_clean[[ll]]
  # sCoor_v <- coordinatesStart_clean[[ll]]
  # eCoor_v <- coordinatesEnd_clean[[ll]]
  # cellID <- cellNames_clean[[ll]]

  nbSNPs <- c()
  indx <- which(chunkL_v < fragLengthToExclude)
  
  if(length(indx) > 0) {
    if(is.element(1, indx)) {
      indx <- indx[-1]
    } else if(is.element(length(chunkL_v), indx)) {
      indx <- indx[-length(chunkL_v)]
    }
    if(length(indx) > 0){
      for(i in indx) {
        #stateSequences_clean[[which(chunksLengths_clean[[4]] < fragLengthToExclude)]]
        x0 <- sCoor_v[i]
        x1 <- eCoor_v[i]
        
        nbSNPs <- c(nbSNPs,(length(which(as.numeric(names(hmm_v)) >= x0 & as.numeric(names(hmm_v)) <= x1))))
      }
    } else {
      nbSNPs <- NA
    }
  } else {
    nbSNPs <- NA
  }
  return(nbSNPs)
}
#-----------

########################################################
##########################################
################################
##################### 
#############
#####
#                        ------------- END OF FUNCTIONS




##########################################
##########################################
##
##
## MAIN FUNCTION
##
##
##########################################
##########################################


wrkDir <- "/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/GT.fullMatrices.Examples"
#input
files <- list.files(path=wrkDir, pattern="*.GT.FORMAT")
openGT <- read.table(paste0(wrkDir, "/", files[15]), header = T)

idFlag <- as.vector(sapply(colnames(openGT)[-c(1:2)], function(t) { unlist(strsplit(t, split = "\\."))[1] } ))
colnames(openGT) <- c(colnames(openGT)[1:2], paste0(idFlag, ".",paste0("SC",1:(ncol(openGT)-2)), "_", paste0("SC",1:(ncol(openGT)-2))))
#output
output_prefix <- sub(".GT.FORMAT", "", files[15])

#getting ID tag 
indTAG <- unlist(strsplit(output_prefix, split = "\\."))[1]
chrID <- unlist(strsplit(output_prefix, split = "\\."))[3]
if(!dir.exists(paste0(wrkDir, "/", indTAG))) {
  print(paste0("Creating directory for individual: ", indTAG))
  dir.create(paste0(wrkDir, "/", indTAG))
  if(!dir.exists(paste0(wrkDir, "/", indTAG, "/", chrID))) {
    dir.create(paste0(wrkDir, "/", indTAG, "/", chrID))
    print(paste0("Creating dir for chrom: ", chrID))
  } else {
    print(paste0("Directory for chrom: ", chrID, " already exists."))
  }
} else {
  print(paste0("Directory: ", indTAG, " already exists"))
  if(!dir.exists(paste0(wrkDir, "/", indTAG, "/", chrID))) {
    dir.create(paste0(wrkDir, "/", indTAG, "/", chrID))
    print(paste0("Creating dir for chrom: ", chrID))
  } else {
    print(paste0("Directory for chrom: ", chrID, " already exists."))
  }
}
## FILTER GT matrix 
##removing cells with low and large number of var observed
nb_var_cell <- unlist(lapply(3:ncol(openGT), function(c) {length(which(openGT[,c] != "."))}))
cells_to_rm <- sort(c(which(nb_var_cell < 50)+2, which(nb_var_cell >= mean(nb_var_cell)+3*sd(nb_var_cell))+2), decreasing = F)
if(length(cells_to_rm) > 0) {
  write.table(matrix(names(openGT)[cells_to_rm],ncol=1), file=paste0(wrkDir,"/", indTAG, "/", chrID, "/", output_prefix, ".excludedCells.txt"), quote = F, row.names = F, col.names = F)
  openGT <- openGT[,-cells_to_rm]
}

##summary stats - cells covering a var
#removes sites with very few cells covering it
nb_cell_var <- unlist(lapply(3:nrow(openGT), function(r) {length(which(openGT[r,3:ncol(openGT)] != "."))}))
vars_to_rm <- sort(which(nb_cell_var <= mean(nb_cell_var)-3*sd(nb_cell_var)), decreasing = F)
if(length(vars_to_rm) > 0) {
  openGT <- openGT[-vars_to_rm,]
}

dim(openGT)
coord_start_end <- c(openGT[1,2],openGT[nrow(openGT),2])
minCoord <- coord_start_end[1]
windowLength <- 150
varToKeep_afterPruning <- c()

while(T) {
  print(minCoord)
  keep <- c()
  if((minCoord+(windowLength-1)) <= coord_start_end[2]) {
    maxCoord <- minCoord+(windowLength-1)
    
  } else {
    maxCoord <- coord_start_end[2]
  }
  v_Pos_tmp <- minCoord:maxCoord
  indexMatching <- intersect(v_Pos_tmp, openGT[,2])

  if(length(indexMatching) > 1) {
    keep <- sample(indexMatching, 1)
  } else if (length(indexMatching) == 1) {
    keep <- indexMatching
  } 
  varToKeep_afterPruning <- c(varToKeep_afterPruning, keep)
  minCoord <- minCoord+windowLength  
  if(minCoord+(windowLength-1) > coord_start_end[2]) {
    break;
  }
}
openGT <- openGT[match(varToKeep_afterPruning, openGT[,2]),]
dim(openGT)
## END of FILTER GT matrix 

##################################
##################################
##
## This part will print SUMMARY STATISTICS across 
## cells and variants
##
##################################
##################################

#summary statistics - GOOD quality cells and sites
pdf(file = paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".summStats.pdf"), height = 4, width = 8)
par(mfrow=c(1,2))
cat(paste("Total_nb_pos", "Total_nb_cells", "Mean_var_per_cell", "Median_var_per_cell", "Mode_var_per_cell", "Lower_range_var_per_cell", "Higher_range_var_per_cell", 
          "Mean_cells_cov_var",  "Median_cells_cov_var", "Mode_cells_cov_var", "Lower_range_cells_cov_var", "Higher_range_cells_cov_var", sep = "\t"),
    file =  paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".summStats.txt"), sep = "\n")
#nb cells and sites
nb_cells <- ncol(openGT)-2
nb_var <- nrow(openGT)
#cells
mean_var_cell <- round(mean(nb_var_cell), digits = 2)
median_var_cell <- round(median(nb_var_cell), digits = 2)
mode_var_cell <- round(density(nb_var_cell)$x[which(density(nb_var_cell)$y == max(density(nb_var_cell)$y))], digits = 2)
range_var_cell <- range(nb_var_cell)
plot(density(nb_var_cell), xlab="Nb of variants covered per cell", main = paste0("Total number of hetSNPs = ", nb_var))

#vars
mean_cells_cov_var <- round(mean(nb_cell_var), digits = 2)
median_cells_cov_var <- round(median(nb_cell_var), digits = 2)
mode_cells_cov_var <- round(density(nb_cell_var)$x[which(density(nb_cell_var)$y == max(density(nb_cell_var)$y))], digits = 2)
range_cells_cov_var <- range(nb_cell_var)
plot(density(nb_cell_var), xlab="Nb of cells covering a var", main = paste0("Total number of cells = ", nb_cells))

cat(paste(nb_var, nb_cells, mean_var_cell, median_var_cell, mode_var_cell, range_var_cell[1], range_var_cell[2], 
          mean_cells_cov_var, median_cells_cov_var, mode_cells_cov_var, range_cells_cov_var[1], 
          range_cells_cov_var[2], sep = "\t"),
    file =  paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".summStats.txt"), sep = "\n", append = T)

dev.off()

## END of SUMMARY STATISTICS

##################################
##################################
##
## This part will compute for each cell the similarity with other cells.
## Similarity here means how many sites in a target cell match 100% of the sites covered by other cells
## or do not match at all other cells
##
##################################
##################################

#declaring variables
nb_cells <- ncol(openGT)-2
nb_var <- nrow(openGT)

targetCell <- 1 #its variable in the loop. it indicates which col to put in 3rd place in the new matrix

while(targetCell <= (ncol(openGT)-2)) {
  
  list_match <- list()
  list_common_var <- list()
  nameCell <- unlist(strsplit(colnames(openGT)[targetCell+2], split = "\\."))[2]
  
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".", nameCell,".sc_comparison.pdf"),  height = 4, width = 8)
  par(mfrow=c(2,2))
  
  comparingCell <- 2 #always fixed
  
  tmpOpenGT <- openGT[,c(1,2,targetCell+2,(3:ncol(openGT))[-targetCell])]
  
  while(comparingCell <= (ncol(tmpOpenGT)-2)) {
    
    #retrieving positions with INFORMATION for target cell
    cell_v_noMissing <- tmpOpenGT[,3][which(tmpOpenGT[,3] != ".")] #target cell will always be in col #3
    names(cell_v_noMissing) <- tmpOpenGT[which(tmpOpenGT[,3] != "."),2]
    nb_var_INFO_target <- length(cell_v_noMissing)
    
    #retrieving positions with INFORMATION for comparing cell
    cell_v_noMissing_compCell <- tmpOpenGT[,2+comparingCell][which(tmpOpenGT[,2+comparingCell] != ".")]
    names(cell_v_noMissing_compCell) <- tmpOpenGT[which(tmpOpenGT[,2+comparingCell] != "."),2]
    
    # overlapping informative sites between the target and the remaining cells.
    overSites <- intersect(names(cell_v_noMissing), names(cell_v_noMissing_compCell))
    list_common_var[[comparingCell]] <- length(overSites)
    if(length(overSites) > 10) {
      match_v <- as.vector(tmpOpenGT[match(overSites,tmpOpenGT[,2]),c(3,comparingCell+2)][,1]) ==  as.vector(tmpOpenGT[match(overSites,tmpOpenGT[,2]),c(3,comparingCell+2)][,2])  
      names(match_v) <- overSites
      list_match[[comparingCell]] <- match_v
    } else {
      list_match[[comparingCell]] <- NA
    }

    comparingCell <- comparingCell+1
  }
  
  #computing similarity, ie number of matching sites
  similarity <- unlist(lapply(list_match[-1], function(comp) {
    sum(comp)/length(comp)
  }))
  hist(similarity, main = nameCell)
  
  #computing the number of informative sites across pairwise combinations of cells
  nb_informative_sites <- unlist(list_common_var)
  hist(nb_informative_sites/nb_var_INFO_target, xlab="% of infor. sites per cell pair", main = nameCell)
  prop_of_similar_cells <- sum(similarity < 0.05 | similarity > 0.95, na.rm = T)/(nb_cells-1)
  #testing whether this cell contains the parental haplotypes
  # if this is the case one should expect this cell to have around 25% complete matching cells and 
  # 25% of zero matching cells, ie 50%
  #this print a list of cells potentially carrying the parental haplotypes.  
  binomTest <- binom.test(sum(similarity < 0.05 | similarity > 0.95, na.rm = T)+1, nb_cells, p = 0.5)
  if(binomTest$p.value > 0.05) {
  
    cat(nameCell, file = paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".", nameCell,".parentalCells.txt"), sep = "\n")
    if(length(which(similarity < 0.05 | similarity > 0.95)) > 0) {
      for(cell in which(similarity < 0.05 | similarity > 0.95)) {
        tmp_cellName <- unlist(strsplit(colnames(tmpOpenGT)[cell+3], split = "\\."))[2]
        cat(paste(tmp_cellName, similarity[cell], sep = "\t"), file = paste0(wrkDir, "/", indTAG, "/", chrID, "/", output_prefix, ".", nameCell,".parentalCells.txt"), sep = "\n", append = T)
      }
    }
  }
  
  #seq_v contains two lists: 1) the number of seq that are different from the target cell and the seq itself.
  seq_v <- lapply(list_match[-1], FUN = retrieve_sequences)
  seqList <- sapply(seq_v, function(seq) {seq$seq})
  seqNb <- sapply(seq_v, function(seq) {seq$nb_of_seq})
  
  new_list_matches_and_errorsNb <- lapply(2:length(list_match), FUN=remove_errors_from_list_match)
  newMatchSeq <- sapply(new_list_matches_and_errorsNb, function(seq) {seq$cleanMatches})
  nbOfErrors <- unlist(sapply(new_list_matches_and_errorsNb, function(seq) {seq$errorNb}))
  mean(nbOfErrors/lengths(list_match[-1]), na.rm = T)
  
  table(lengths(seqList))
  hist(lengths(seqList), breaks = 0:max(lengths(seqList)), xlab = "Recombining length (nbSNPs)", main = nameCell)
  hist(seqNb, breaks = -0.5:(max(seqNb, na.rm = T)+0.5), xlab = "Nb of recombining segments", main = nameCell)
  dev.off()
 
  targetCell <- targetCell+1 
}

#####################--------------------------

##################################
##################################
##
## This part will take the cells identified to represent 
## the parental chromosomes and will compute
## the proportion of errors
##
##################################
##################################
require(filesstrings)
pathToIndDir  <- paste0(wrkDir, "/", indTAG, "/", chrID)

if(!dir.exists(pathToIndDir)) {

  dir.create(pathToIndDir, showWarnings = F)
  listOfFilesInd <- list.files(path=paste0(wrkDir, "/"), 
                               pattern = paste0(output_prefix, "*"))
  fileToKeep <- list.files(path=paste0(wrkDir, "/"), pattern=paste0(indID, ".*.GT.FORMAT"))
  listOfFilesIndToCopy <- setdiff(listOfFilesInd, fileToKeep)
  lapply(listOfFilesIndToCopy, function(x) { file.move(paste0(wrkDir, "/", x), pathToIndDir) }) 
  
}

parentalCells <- list.files(path=pathToIndDir, pattern="*.parentalCells.txt")

if(length(parentalCells) > 10) { #the remaining code will be run ONLY and ONLY if there are at least 10 parental cells.

  listParentCells <- as.character(sapply(parentalCells, function(x) { unlist(strsplit(x, split = "\\."))[4] }))
  indxParentCells <- as.numeric(sapply(listParentCells, function(x) { grep(x, colnames(openGT))}))
  
  #here the openGT is already without bad quality cells and variants. 
  newOpenGT <- openGT[,c(1,2,indxParentCells)]
  similarity_m <- matrix(rep(NA, length(indxParentCells)*length(indxParentCells)), ncol = length(indxParentCells))
  colnames(similarity_m) <- listParentCells
  COUNT_CELLS <- 1
  for(file in parentalCells) {
    
    #file <- parentalCells[1]
    readFile <- read.table(paste0(pathToIndDir, "/", file), skip = 1)
    interCells <- intersect(readFile$V1, colnames(similarity_m))
    similarity_m[COUNT_CELLS,match(interCells, colnames(similarity_m))] <- readFile$V2[match(interCells,readFile$V1)]
    COUNT_CELLS <- COUNT_CELLS+1
  }
  df_similarity <- similarity_m[which(is.na(similarity_m))]  <- 0.5
  df_similarity <- scale(similarity_m)
  kmeans_similarity <- kmeans(df_similarity, centers = 2)
  parentalONE <- colnames(df_similarity)[kmeans_similarity$cluster == 2] 
  parentalTWO <- colnames(df_similarity)[kmeans_similarity$cluster == 1] 
  
  newOpenGT[,match(parentalONE, listParentCells)+2]
  cat(output_prefix, sep="\n", file=paste0(pathToIndDir, "/", output_prefix, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"))
  cat(output_prefix, sep="\n", file=paste0(pathToIndDir, "/", output_prefix, ".errorRatePerSNP.txt"))
  alleleStates <- c("0","1")
  haplotypes <- c()
  for(i in 1:nrow(newOpenGT)) {
    
    #i <- 1
    #parental ONE have equal nb of cells carrying each of the alleles 
    print(i)
    error <- NA
    #getting the table of alleles contained by each of the parentals (alleles present in each set and how many cells carrying those same alleles)
    #in a perfect world each of the parental sets should ALL have the same allele (often this is not the case)
    pONE <- sum(table(newOpenGT[i,match(parentalONE, listParentCells)+2][newOpenGT[i,match(parentalONE, listParentCells)+2] != "."] ))
    pTWO <- sum(table(newOpenGT[i,match(parentalTWO, listParentCells)+2][newOpenGT[i,match(parentalTWO, listParentCells)+2] != "."] ))
  
    # print(table(newOpenGT[i,match(parentalONE, listParentCells)+2][newOpenGT[i,match(parentalONE, listParentCells)+2] != "."] ))
    # print(table(newOpenGT[i,match(parentalTWO, listParentCells)+2][newOpenGT[i,match(parentalTWO, listParentCells)+2] != "."] ))
    
    if(sum(pONE,pTWO) < 3) { #less than 3 cells cover the site ---> umbiguous state
      alleleONE <- NA
      alleleTWO <- NA
      #print(paste0(alleleONE, alleleTWO))
    } else { #at least three cells cover a site
      if(pONE > 0 & pTWO == 0) {  #there are at least three cells covering one site in group ONE of cells, we can have one or two alleles
  
        summaryONE <- table(newOpenGT[i,match(parentalONE, listParentCells)+2][newOpenGT[i,match(parentalONE, listParentCells)+2] != "."] )
        maxAlleleONE <- names(summaryONE)[summaryONE == max(summaryONE)]
  
        if(length(summaryONE) == 2){     #if we have two alleles in the set ONE of cells
          
          ratio_twoAlleles <- max(summaryONE/sum(summaryONE))
          
          if(ratio_twoAlleles >= (4/5)) { #if one of the two alleles is present more than 4/5 times than the other allele
            
            alleleONE <- maxAlleleONE
            alleleTWO <- setdiff(alleleStates, maxAlleleONE)
              
          } else { #the ratio of the diff in freq is less then 4/5
            
            alleleONE <- NA
            alleleTWO <- NA
          }
  
        } else if(length(summaryONE) == 1) { #if we have ONLY ONE allele in the set ONE of cells
  
          if(summaryONE[summaryONE == max(summaryONE)] >= 5) { #the allele is seen at least 5 times
            
            alleleONE <- maxAlleleONE
            alleleTWO <- setdiff(alleleStates, maxAlleleONE)
            
          } else { #the allele is seen less than 5 times
            
            alleleONE <- NA
            alleleTWO <- NA
          }
        }
      } else if(pONE == 0 & pTWO > 0) { #there are at least five cells covering one site in group TWO of cells, we can have one or two alleles
  
        summaryTWO <- table(newOpenGT[i,match(parentalTWO, listParentCells)+2][newOpenGT[i,match(parentalTWO, listParentCells)+2] != "."] )
        maxAlleleTWO <- names(summaryTWO)[summaryTWO == max(summaryTWO)]
        
        if(length(summaryTWO) == 2){ #if we have two alleles in the set TWO of cells
          
          ratio_twoAlleles <- max(summaryTWO/sum(summaryTWO))
          
          if(ratio_twoAlleles >= (4/5)) { #if one of the two alleles is present more than 4/5 times than the other allele
            
            alleleTWO <- maxAlleleTWO
            alleleONE <- setdiff(alleleStates, maxAlleleTWO)
            
          } else { #the ratio of the diff in freq is less then 4/5
            
            alleleONE <- NA
            alleleTWO <- NA
          }
          
        } else if(length(summaryTWO) == 1) { #if we have ONLY ONE allele in the set TWO of cells
          
          if(summaryTWO[summaryTWO == max(summaryTWO)] >= 5) { #the allele is seen at least 5 times
            
            alleleTWO <- maxAlleleTWO
            alleleONE <- setdiff(alleleStates, maxAlleleTWO)
            
          } else { #the allele is seen less than 5 times
            
            alleleONE <- NA
            alleleTWO <- NA
          }
        }
        
      } else if(pONE > 0 & pTWO > 0) { #there are at least 5 cells in both set ONE and TWO covering a site 
  
        summaryONE <- table(newOpenGT[i,match(parentalONE, listParentCells)+2][newOpenGT[i,match(parentalONE, listParentCells)+2] != "."] )
        maxAlleleONE <- names(summaryONE)[summaryONE == max(summaryONE)]
        
        summaryTWO <- table(newOpenGT[i,match(parentalTWO, listParentCells)+2][newOpenGT[i,match(parentalTWO, listParentCells)+2] != "."] )
        maxAlleleTWO <- names(summaryTWO)[summaryTWO == max(summaryTWO)]
        
        if(length(summaryONE) == 1 & length(summaryTWO) == 1) { #both groups have ONLY one allele among the cells
          
          if(maxAlleleONE == maxAlleleTWO) {
            
            alleleONE <- NA
            alleleTWO <- NA
            
          } else { #the alleles present in the two sets of cells are different (regardeless of their ind. frequency)
           
            alleleONE <- maxAlleleONE
            alleleTWO <- maxAlleleTWO
            
          }
  
        } else if(length(summaryONE) == 1 & length(summaryTWO) == 2) { #both groups have TWO alleles among the cells ---> umbiguous
  
          ratio_twoAlleles <- max(summaryTWO/sum(summaryTWO))
  
          if(length(maxAlleleTWO) == 2) {
            alleleONE <- NA
            alleleTWO <- NA
          
          } else {
            if(ratio_twoAlleles >= (4/5) & maxAlleleONE != maxAlleleTWO) { #if one of the two alleles is present more than 4/5 times than the other allele
              
              alleleTWO <- maxAlleleTWO
              alleleONE <- maxAlleleONE
              #error allele
              error <- summaryTWO[which(names(summaryTWO) != maxAlleleTWO)]/pTWO
              
            } else { #the ratio of the diff in freq is less then 4/5
              
              alleleONE <- NA
              alleleTWO <- NA
            }
          }
        } else if(length(summaryONE) == 2 & length(summaryTWO) == 1) { #set ONE has two alleles and set TWO has one 
          
          ratio_twoAlleles <- max(summaryONE/sum(summaryONE))
          
          if(length(maxAlleleONE) == 2) { #alleles in set ONE have the same frequency across cells ----- > ambiguous
            
            alleleONE <- NA
            alleleTWO <- NA
            
          } else {
            if(ratio_twoAlleles >= (4/5) & maxAlleleONE != maxAlleleTWO) { #if one of the two alleles is present more than 4/5 times than the other allele
              
              alleleONE <- maxAlleleONE
              alleleTWO <- maxAlleleTWO
              #error allele
              error <- summaryONE[which(names(summaryONE) != maxAlleleONE)]/pONE
              
            } else { #the ratio of the diff in freq is less then 4/5
              
              alleleONE <- NA
              alleleTWO <- NA
            }          
          }
        } else if(length(summaryONE) == 2 & length(summaryTWO) == 2) { #set ONE and set TWO have two alleles
  
          if(length(maxAlleleTWO) == 1 & length(maxAlleleONE) == 1) {
            
            ratio_twoAllelesONE <- max(summaryONE/sum(summaryONE))      
            ratio_twoAllelesTWO <- max(summaryTWO/sum(summaryTWO))
  
            if(ratio_twoAllelesONE >= (4/5) & ratio_twoAllelesTWO >= (4/5) & maxAlleleONE != maxAlleleTWO) { #in both sets the most freq allele is at a freq > 4/5 
              # also the most freq allele of set ONE is NOT the most freq allele of set TWO
              alleleONE <- maxAlleleONE
              alleleTWO <- maxAlleleTWO
              #error allele
              errorONE <- summaryONE[which(names(summaryONE) != maxAlleleONE)]/pONE
              errorTWO <- summaryTWO[which(names(summaryTWO) != maxAlleleTWO)]/pTWO
              error <- errorONE+errorTWO
  
            } else {
              
              alleleONE <- NA
              alleleTWO <- NA
              
            } 
          } else { #when either in set ONE or in set TWO the two existing alleles are at the same freq
  
            alleleONE <- NA
            alleleTWO <- NA
            
          }
        } #end of the case where there are two alleles in both sets of cells
      } #end of the cases where there is INFO in both sets of cells
    } #end of the cases where there is at least 5 cells across the whole set of cells covering a site.  
    print(c(newOpenGT[i,2],alleleONE, alleleTWO))
    print(c(newOpenGT[i,2], error))
    cat(paste(newOpenGT[i,2],alleleONE, alleleTWO, sep = "\t"), sep="\n", file=paste0(pathToIndDir, "/", output_prefix, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"),
        append = T)
    cat(paste(newOpenGT[i,2],error, sep = "\t"), sep="\n", file=paste0(pathToIndDir, "/", output_prefix, ".errorRatePerSNP.txt"), append = T)
    # errorAlleleONE <-  
    # errorAlleleTWO <- 
      
    #haplotypes <- rbind(haplotypes, c(newOpenGT[i,2],alleleONE, alleleTWO))
  }
  
  
  hap <- read.table(file=paste0(pathToIndDir, "/", output_prefix, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), skip = 1, header = F, row.names = 1)
  hap <- hap[which(!is.na(hap$V2)),]
  inverted_hap <- t(hap)
  write.table(inverted_hap, file = paste0(pathToIndDir, "/", indTAG, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  
  ##################################
  ##################################
  ##
  ## This plots the raw bp by comparing the 
  ## inferred parental haplotypes with every single
  ## cell in windows of 10 SNPs.
  ##
  ##
  ##################################
  ##################################
  
  #opening haplotype file
  inverted_hap <- read.table(paste0(pathToIndDir, "/", indTAG, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), header = T)
  hapFile_posNames <- scan(paste0(pathToIndDir, "/", indTAG, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), nlines=1, what=numeric())
  openHap <- inverted_hap
  colnames(openHap) <- hapFile_posNames
  
  outfileName <- paste0(pathToIndDir, "/SC_phasedHapHMap_", indTAG,"_",chrID, "_QCvar25_NbCells10_NbLinks5.pdf") 
  
  #here openGT is already filtered out
  subSet_SC <- as.matrix(openGT[,3:ncol(openGT)])
  subSet_SC[which(subSet_SC == ".")] <- NA
  subSet_SC <- as.data.frame(subSet_SC)
  subSet_SC <- as.matrix(subSet_SC[match(colnames(openHap), openGT[,2]),])
  subVarNames <- openGT[,2][match(colnames(openHap), openGT[,2])]
  
  #S23_chr9_QCvar25_NbCells10_NbLinks5.pdf
  pdf(file=outfileName, height = 11, width = 6)
  par(mfrow=c(5,1))
  cellNb <- 1
  
  while (cellNb <= ncol(subSet_SC)) {
    
    #cellNb <- 1
    index_to_compare <- which(!is.na(subSet_SC[,cellNb]))
    window <- 1
    
    
    hapOne <- c()
    hapTwo <- c()
    listWindows <- list()
    
    while (window+9 <= length(index_to_compare)) {
      
      vectPerCell <- subSet_SC[index_to_compare[window:(window+9)],cellNb]
      namesWindow <- subVarNames[index_to_compare[window:(window+9)]]
      names(vectPerCell) <- namesWindow
      
      subHapList <- openHap[,match(namesWindow, as.numeric(colnames(openHap)))]
      
      hapOne[window] <- sum(unlist(lapply(1:length(vectPerCell), function(x) { vectPerCell[x] == subHapList[1,x]  })))/length(vectPerCell)
      hapTwo[window] <- sum(unlist(lapply(1:length(vectPerCell), function(x) { vectPerCell[x] == subHapList[2,x]  })))/length(vectPerCell)
      
      listWindows[[window]] <- vectPerCell
      window <- window+1
      
    }
    plot(hapOne, type="l", col="blue", ylim=c(0,1))
    lines(hapTwo, type="l", col="red")
    # windowsAcrossCells[[cellNb]] <- listWindows
    # propAcrossCells[[cellNb]] <- rbind(hapOne,hapTwo) 
    cellNb <- cellNb+1
  }
  dev.off()
  ####----------------
  ##----------
  #-----
  
    ##################################
    ##################################
    ##
    ## This part will take the every cell of the filtered 
    ## geno matrix and will apply a hmm to infer the
    ## hidden states. 
    ##
    ##################################
    ##################################
    library(parallel)
    namesParHaps <- scan(paste0(pathToIndDir, "/", indTAG, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), nlines = 1)
    openParentalHaps <- read.table(paste0(pathToIndDir, "/", indTAG, ".finalHapNeighborHaps_HQ25_NbCells10_NbLinks5.txt"), header = T)
    colnames(openParentalHaps) <- namesParHaps
    ncol(openParentalHaps)
    
    errorRates <- read.table(paste0(pathToIndDir, "/", indTAG, ".SingleCellsHetSNPs.chr20.errorRatePerSNP.txt"), header = F, skip = 1)
    nrow(errorRates)
    
    recombRate_perSite <- 1.25*10^-8
    states <- c("P", "M")
    symbols <- c("0", "1")
    startProb <- 0.5
    
    dim(openGT)
    no_cores <- detectCores()-1
    cl <- makeCluster(no_cores)
    clusterExport(cl, c("openGT", "inferring_bp_hmm", "openParentalHaps", "errorRates", "recombRate_perSite", "states", "symbols", "startProb"))
    hmm_results <- parLapply(cl,1:(ncol(openGT)-2), function(x) { inferring_bp_hmm(x) })
    #coco <- lapply(1:(ncol(openGT)-2), function(x) { inferring_bp_hmm(x) })
    stopCluster(cl)
    
    
    cellNames_v <- sapply(hmm_results, function(cell) {cell$cellName})
    cell_hmm_seq_m <- t(do.call(rbind, lapply(hmm_results, function(cell) {cell$inferredSeq})))
    colnames(cell_hmm_seq_m) <- cellNames_v
    cell_hmm_StateSeq_m <- lapply(hmm_results, function(cell) {cell$inferredStates})
    #cell_hmm_StateSeq_m <- t(do.call(rbind, lapply(hmm_results, function(cell) {cell$inferredStates})))
    #colnames(cell_hmm_StateSeq_m) <- cellNames_v
    
    subSet_SC <- cell_hmm_seq_m
    subVarNames <- rownames(cell_hmm_seq_m)
    outfileName <- paste0(pathToIndDir, "/SC_phasedHapHMap_", indTAG,"_",chrID, "_QCvar25_NbCells10_NbLinks5_postHMM.pdf") 
      
    subSet_SC[subSet_SC == "."] <- NA
    #S23_chr9_QCvar25_NbCells10_NbLinks5.pdf
    pdf(file=outfileName, height = 11, width = 6)
    par(mfrow=c(5,1))
    cellNb <- 1
      
    while (cellNb <= ncol(subSet_SC)) {
      
      #cellNb <- 1
      index_to_compare <- which(!is.na(subSet_SC[,cellNb]))
      window <- 1
      
      
      hapOne <- c()
      hapTwo <- c()
      listWindows <- list()
      
      while (window+9 <= length(index_to_compare)) {
        
        vectPerCell <- subSet_SC[index_to_compare[window:(window+9)],cellNb]
        namesWindow <- subVarNames[index_to_compare[window:(window+9)]]
        names(vectPerCell) <- namesWindow
        
        subHapList <- openHap[,match(namesWindow, as.numeric(colnames(openHap)))]
        
        hapOne[window] <- sum(unlist(lapply(1:length(vectPerCell), function(x) { vectPerCell[x] == subHapList[1,x]  })))/length(vectPerCell)
        hapTwo[window] <- sum(unlist(lapply(1:length(vectPerCell), function(x) { vectPerCell[x] == subHapList[2,x]  })))/length(vectPerCell)
        
        listWindows[[window]] <- vectPerCell
        window <- window+1
        
      }
      plot(hapOne, type="l", col="blue", ylim=c(0,1))
      lines(hapTwo, type="l", col="red")
      # windowsAcrossCells[[cellNb]] <- listWindows
      # propAcrossCells[[cellNb]] <- rbind(hapOne,hapTwo) 
      cellNb <- cellNb+1
    }
    dev.off()
  
  ###############################################
  ##
  ## plot inferred parental blocks 
  ##
  ###############################################
  #####################                                                                                                                                      
  ##
  ## instanciate maternal and paternal cells                                                                                                                                  
  ##
  #####################  
  
  
  finalInference <- cell_hmm_StateSeq_m
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/InferredHapsHMap_", indTAG, "_", chrID,"_QCvar25_NbCells10_NbLinks5.pdf"), height = 9, width = 6)
  
  featuresCells <- lapply(1:length(finalInference), function(cellNb) { names(table(finalInference[[cellNb]]))  })
  cellsWithUnknown <- lapply(1:length(finalInference), function(cellNb) { is.element("U",featuresCells[[cellNb]]) })
  cellsWithOUTUnknown <- lapply(1:length(finalInference), function(cellNb) { !is.element("U",featuresCells[[cellNb]]) })
  
  indexCellsUnknown <- which(cellsWithUnknown == T)
  indexCellsWithOUTUnknown <- which(cellsWithOUTUnknown == T)
  
  onlyMcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("P",featuresCells[[cellNb]]) })
  onlyPcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("M",featuresCells[[cellNb]]) })
  
  indexOnlyPCells <- which(onlyPcells == T)
  indexOnlyMCells <- which(onlyMcells == T)
  
  PandMCells <- lapply(1:length(finalInference), function(cellNb) { is.element("M",featuresCells[[cellNb]]) & is.element("P",featuresCells[[cellNb]]) })
  
  PStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "P")[1] })
  MStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "M")[1] })
  
  doesItStartWithP <- lapply(1:length(finalInference), function(cellNb) { PStartingChr[[cellNb]] < MStartingChr[[cellNb]]  })
  doesItStartWithM <- lapply(1:length(finalInference), function(cellNb) { MStartingChr[[cellNb]] < PStartingChr[[cellNb]]  })
  
  PstrechesSize <- lapply(1:length(finalInference), function(cellNb) { abs(PStartingChr[[cellNb]]-MStartingChr[[cellNb]]) })
  
  startP_decreasing <- which(PandMCells == T & doesItStartWithP == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithP == T))]), decreasing = T)]
  startM_decreasing <- which(PandMCells == T & doesItStartWithM == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithM == T))]), decreasing = T)]
  
  
  finalPaternalCells <- c(indexOnlyPCells,startP_decreasing)
  finalMaternalCells <- c(indexOnlyMCells,startM_decreasing)
  
  #plot
  plot("", xlim=c(0,length(finalInference[[1]])), ylim = c(0,length(finalInference)), xlab = "hetSNPs", ylab = "Cells")
  nbPolygons <- list()
  
  countSC <- 1
  for (cellNb in c(finalPaternalCells, finalMaternalCells)) {
    
    #cellNb <- 159
    v <- finalInference[[cellNb]]
    
    mCoord <- c()
    pCoord <- c()
    uCoord <- c()
    
    state <- "I"
    for (i in 1:length(v)) {
      
      
      if (v[i] != "." & v[i] == "M" & v[i] != state) {
        
        state <- "M"
        mCoord <- c(mCoord, i)
      } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
        state <- "P"
        pCoord <- c(pCoord, i)
        
      } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
        state <- "U"
        uCoord <- c(uCoord, i)
        
      } 
      
    }
    lastCoord <- as.numeric(which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))])
    x_coord <- sort(c(mCoord, pCoord, uCoord, lastCoord))
    nbPolygons[[cellNb]] <- length(x_coord)-1
    
    
    for (polNb in 1:nbPolygons[[cellNb]]) {
      
      polygon(x=c(x_coord[polNb],x_coord[polNb],x_coord[polNb+1],x_coord[polNb+1]), y=c((countSC-1),countSC,countSC,(countSC-1)), col=colPolygon(x_coord[polNb]), border = NA)  
    }
    countSC <- countSC+1
  }
  dev.off()
  
  ###############################################################
  ##
  ##
  ## ------------ IDENTIFYING PARENTAL CHUNKS -------------------
  ##
  ##
  ###############################################################
  
  list_chr_chunks <- lapply(1:length(cell_hmm_StateSeq_m), function(z) {compute_nb_chromosome_chunks(z)})
  chunksLengths <- lapply(list_chr_chunks, function(z) { z$chunkLength})
  stateSequences <- lapply(list_chr_chunks, function(z) { z$chunksIdentity})
  coordinates_start <- lapply(list_chr_chunks, function(z) { z$sCoordinate})
  coordinates_end <- lapply(list_chr_chunks, function(z) { z$eCoordinate})
  
  list_hmm_clean <- cell_hmm_StateSeq_m[-which(lengths(chunksLengths) > 8)]
  cellNames_clean <- cellNames_v[-which(lengths(chunksLengths) > 8)]
  chunksLengths_clean <- chunksLengths[-which(lengths(chunksLengths) > 8)]
  stateSequences_clean <- stateSequences[-which(lengths(chunksLengths) > 8)]
  coordinatesStart_clean <- coordinates_start[-which(lengths(chunksLengths) > 8)]
  coordinatesEnd_clean <- coordinates_end[-which(lengths(chunksLengths) > 8)]
  
  # length(list_hmm_clean)
  # length(cellNames_clean)
  # length(chunksLengths_clean)
  # length(coordinatesStart_clean)
  # length(coordinatesEnd_clean)
  
  
  #writing to a file cells with more than 8 chunks
  write.table(matrix(cellNames_v[which(lengths(chunksLengths) > 8)]
                     , ncol=1), file = paste0(wrkDir, "/", indTAG, "/", chrID, "/ExcludedCells_highNbofBP.txt"), quote = F, row.names = F, col.names = F)
  
  #write a table with state sequences per cell
  stateSeqHMM_clean <- do.call(rbind, list_hmm_clean)
  rownames(stateSeqHMM_clean)  <- cellNames_clean
  write.table(t(stateSeqHMM_clean), file = paste0(wrkDir, "/", indTAG, "/", chrID, "/",indTAG, ".", chrID, ".stateSeqHMM_clean.txt"), quote = F, 
              row.names = T, col.names = T, sep = "\t")
  
  #################################################
  ## plotting parental chunks without cells with more than 8 chunks
  finalInference <- cell_hmm_StateSeq_m[-which(lengths(chunksLengths) > 8)]
  
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/InferredHapsHMap_", indTAG, "_", chrID,"_QCvar25_NbCells10_NbLinks5_filtered.pdf"), height = 8, width = 6)
  
  featuresCells <- lapply(1:length(finalInference), function(cellNb) { names(table(finalInference[[cellNb]]))  })
  cellsWithUnknown <- lapply(1:length(finalInference), function(cellNb) { is.element("U",featuresCells[[cellNb]]) })
  cellsWithOUTUnknown <- lapply(1:length(finalInference), function(cellNb) { !is.element("U",featuresCells[[cellNb]]) })
  
  indexCellsUnknown <- which(cellsWithUnknown == T)
  indexCellsWithOUTUnknown <- which(cellsWithOUTUnknown == T)
  
  onlyMcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("P",featuresCells[[cellNb]]) })
  onlyPcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("M",featuresCells[[cellNb]]) })
  
  indexOnlyPCells <- which(onlyPcells == T)
  indexOnlyMCells <- which(onlyMcells == T)
  
  PandMCells <- lapply(1:length(finalInference), function(cellNb) { is.element("M",featuresCells[[cellNb]]) & is.element("P",featuresCells[[cellNb]]) })
  
  PStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "P")[1] })
  MStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "M")[1] })
  
  doesItStartWithP <- lapply(1:length(finalInference), function(cellNb) { PStartingChr[[cellNb]] < MStartingChr[[cellNb]]  })
  doesItStartWithM <- lapply(1:length(finalInference), function(cellNb) { MStartingChr[[cellNb]] < PStartingChr[[cellNb]]  })
  
  PstrechesSize <- lapply(1:length(finalInference), function(cellNb) { abs(PStartingChr[[cellNb]]-MStartingChr[[cellNb]]) })
  
  startP_decreasing <- which(PandMCells == T & doesItStartWithP == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithP == T))]), decreasing = T)]
  startM_decreasing <- which(PandMCells == T & doesItStartWithM == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithM == T))]), decreasing = T)]
  
  
  finalPaternalCells <- c(indexOnlyPCells,startP_decreasing)
  finalMaternalCells <- c(indexOnlyMCells,startM_decreasing)
  
  #plot
  plot("", xlim=c(0,length(finalInference[[1]])), ylim = c(0,length(finalInference)), xlab = "hetSNPs", ylab = "Cells")
  nbPolygons <- list()
  
  countSC <- 1
  for (cellNb in c(finalPaternalCells, finalMaternalCells)) {
    
    #cellNb <- 159
    v <- finalInference[[cellNb]]
    
    mCoord <- c()
    pCoord <- c()
    uCoord <- c()
    
    state <- "I"
    for (i in 1:length(v)) {
      
      
      if (v[i] != "." & v[i] == "M" & v[i] != state) {
        
        state <- "M"
        mCoord <- c(mCoord, i)
      } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
        state <- "P"
        pCoord <- c(pCoord, i)
        
      } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
        state <- "U"
        uCoord <- c(uCoord, i)
        
      } 
      
    }
    lastCoord <- as.numeric(which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))])
    x_coord <- sort(c(mCoord, pCoord, uCoord, lastCoord))
    nbPolygons[[cellNb]] <- length(x_coord)-1
    
    
    for (polNb in 1:nbPolygons[[cellNb]]) {
      
      polygon(x=c(x_coord[polNb],x_coord[polNb],x_coord[polNb+1],x_coord[polNb+1]), y=c((countSC-1),countSC,countSC,(countSC-1)), col=colPolygon(x_coord[polNb]), border = NA)  
    }
    countSC <- countSC+1
  }
  dev.off()
  
  #######################
  ##
  ## END OF PLOT CHROM CHUNKS
  ##
  #######################
  
  
  ####################################
  ##
  ##
  ## Cleaning the fragments shorter than a given size
  ##
  ##
  ####################################
  fragLengthToExclude <- 1000000
  
  new_hmm_l <- list()
  snps_per_chunk  <- list()
  
  new_hmm_l <- lapply(1:length(chunksLengths_clean), function(z) {remove_short_chunks(chunksLengths_clean[[z]], list_hmm_clean[[z]], coordinatesStart_clean[[z]], coordinatesEnd_clean[[z]])})
  length(new_hmm_l)
  
  snps_per_chunk <- unlist(lapply(1:length(chunksLengths_clean), function(z) {counting_SNPs_per_chunk(chunksLengths_clean[[z]], list_hmm_clean[[z]], coordinatesStart_clean[[z]], coordinatesEnd_clean[[z]], cellNames_clean[[z]])}))
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/InferredHapsHMap_", indTAG, "_", chrID,"_nbSNPschunk_1MbRemovedChunks.pdf"), height = 5, width = 5)
  hist(snps_per_chunk, breaks = 1:max(snps_per_chunk, na.rm = T))
  dev.off()
  finalInference <- new_hmm_l
  
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/InferredHapsHMap_", indTAG, "_", chrID,"_QCvar25_NbCells10_NbLinks5_filtered_1MbRemoved_new.pdf"), height = 8, width = 6)
  
  featuresCells <- lapply(1:length(finalInference), function(cellNb) { names(table(finalInference[[cellNb]]))  })
  cellsWithUnknown <- lapply(1:length(finalInference), function(cellNb) { is.element("U",featuresCells[[cellNb]]) })
  cellsWithOUTUnknown <- lapply(1:length(finalInference), function(cellNb) { !is.element("U",featuresCells[[cellNb]]) })
  
  indexCellsUnknown <- which(cellsWithUnknown == T)
  indexCellsWithOUTUnknown <- which(cellsWithOUTUnknown == T)
  
  onlyMcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("P",featuresCells[[cellNb]]) })
  onlyPcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("M",featuresCells[[cellNb]]) })
  
  indexOnlyPCells <- which(onlyPcells == T)
  indexOnlyMCells <- which(onlyMcells == T)
  
  PandMCells <- lapply(1:length(finalInference), function(cellNb) { is.element("M",featuresCells[[cellNb]]) & is.element("P",featuresCells[[cellNb]]) })
  
  PStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "P")[1] })
  MStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "M")[1] })
  
  doesItStartWithP <- lapply(1:length(finalInference), function(cellNb) { PStartingChr[[cellNb]] < MStartingChr[[cellNb]]  })
  doesItStartWithM <- lapply(1:length(finalInference), function(cellNb) { MStartingChr[[cellNb]] < PStartingChr[[cellNb]]  })
  
  PstrechesSize <- lapply(1:length(finalInference), function(cellNb) { abs(PStartingChr[[cellNb]]-MStartingChr[[cellNb]]) })
  
  startP_decreasing <- which(PandMCells == T & doesItStartWithP == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithP == T))]), decreasing = T)]
  startM_decreasing <- which(PandMCells == T & doesItStartWithM == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithM == T))]), decreasing = T)]
  
  
  finalPaternalCells <- c(indexOnlyPCells,startP_decreasing)
  finalMaternalCells <- c(indexOnlyMCells,startM_decreasing)
  
  #plot
  plot("", xlim=c(0,length(finalInference[[1]])), ylim = c(0,length(finalInference)), xlab = "hetSNPs", ylab = "Cells")
  nbPolygons <- list()
  
  countSC <- 1
  for (cellNb in c(finalPaternalCells, finalMaternalCells)) {
    
    #cellNb <- 159
    v <- finalInference[[cellNb]]
    
    mCoord <- c()
    pCoord <- c()
    uCoord <- c()
    
    state <- "I"
    for (i in 1:length(v)) {
      
      
      if (v[i] != "." & v[i] == "M" & v[i] != state) {
        
        state <- "M"
        mCoord <- c(mCoord, i)
      } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
        state <- "P"
        pCoord <- c(pCoord, i)
        
      } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
        state <- "U"
        uCoord <- c(uCoord, i)
        
      } 
      
    }
    lastCoord <- as.numeric(which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))])
    x_coord <- sort(c(mCoord, pCoord, uCoord, lastCoord))
    nbPolygons[[cellNb]] <- length(x_coord)-1
    
    
    for (polNb in 1:nbPolygons[[cellNb]]) {
      
      polygon(x=c(x_coord[polNb],x_coord[polNb],x_coord[polNb+1],x_coord[polNb+1]), y=c((countSC-1),countSC,countSC,(countSC-1)), col=colPolygon(x_coord[polNb]), border = NA)  
    }
    countSC <- countSC+1
  }
  dev.off()
  
  #######################
  ##
  ## END OF PLOT CHROM CHUNKS
  ##
  #######################
  
  fragLengthToExclude <- 5000000
  
  new_hmm_l <- lapply(1:length(chunksLengths_clean), function(z) {remove_short_chunks(chunksLengths_clean[[z]], list_hmm_clean[[z]], coordinatesStart_clean[[z]], coordinatesEnd_clean[[z]])})
  length(new_hmm_l)
  
  finalInference <- new_hmm_l
  
  pdf(file=paste0(wrkDir, "/", indTAG, "/", chrID, "/InferredHapsHMap_", indTAG, "_", chrID,"_QCvar25_NbCells10_NbLinks5_filtered_5MbRemoved.pdf"), height = 8, width = 6)
  
  featuresCells <- lapply(1:length(finalInference), function(cellNb) { names(table(finalInference[[cellNb]]))  })
  cellsWithUnknown <- lapply(1:length(finalInference), function(cellNb) { is.element("U",featuresCells[[cellNb]]) })
  cellsWithOUTUnknown <- lapply(1:length(finalInference), function(cellNb) { !is.element("U",featuresCells[[cellNb]]) })
  
  indexCellsUnknown <- which(cellsWithUnknown == T)
  indexCellsWithOUTUnknown <- which(cellsWithOUTUnknown == T)
  
  onlyMcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("P",featuresCells[[cellNb]]) })
  onlyPcells <- lapply(1:length(finalInference), function(cellNb) { !is.element("M",featuresCells[[cellNb]]) })
  
  indexOnlyPCells <- which(onlyPcells == T)
  indexOnlyMCells <- which(onlyMcells == T)
  
  PandMCells <- lapply(1:length(finalInference), function(cellNb) { is.element("M",featuresCells[[cellNb]]) & is.element("P",featuresCells[[cellNb]]) })
  
  PStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "P")[1] })
  MStartingChr <- lapply(1:length(finalInference), function(cellNb) {  which(finalInference[[cellNb]] == "M")[1] })
  
  doesItStartWithP <- lapply(1:length(finalInference), function(cellNb) { PStartingChr[[cellNb]] < MStartingChr[[cellNb]]  })
  doesItStartWithM <- lapply(1:length(finalInference), function(cellNb) { MStartingChr[[cellNb]] < PStartingChr[[cellNb]]  })
  
  PstrechesSize <- lapply(1:length(finalInference), function(cellNb) { abs(PStartingChr[[cellNb]]-MStartingChr[[cellNb]]) })
  
  startP_decreasing <- which(PandMCells == T & doesItStartWithP == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithP == T))]), decreasing = T)]
  startM_decreasing <- which(PandMCells == T & doesItStartWithM == T)[order(unlist(PstrechesSize[c(which(PandMCells == T & doesItStartWithM == T))]), decreasing = T)]
  
  
  finalPaternalCells <- c(indexOnlyPCells,startP_decreasing)
  finalMaternalCells <- c(indexOnlyMCells,startM_decreasing)
  
  #plot
  plot("", xlim=c(0,length(finalInference[[1]])), ylim = c(0,length(finalInference)), xlab = "hetSNPs", ylab = "Cells")
  nbPolygons <- list()
  
  countSC <- 1
  for (cellNb in c(finalPaternalCells, finalMaternalCells)) { #length(new_hmm_l)
    
    #print(cellNb)
    #cellNb <- 159
    v <- finalInference[[cellNb]]
    
    mCoord <- c()
    pCoord <- c()
    uCoord <- c()
    
    state <- "I"
    for (i in 1:length(v)) {
      
      
      if (v[i] != "." & v[i] == "M" & v[i] != state) {
        
        state <- "M"
        mCoord <- c(mCoord, i)
      } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
        state <- "P"
        pCoord <- c(pCoord, i)
        
      } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
        state <- "U"
        uCoord <- c(uCoord, i)
        
      } 
      
    }
    lastCoord <- as.numeric(which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))])
    x_coord <- sort(c(mCoord, pCoord, uCoord, lastCoord))
    nbPolygons[[cellNb]] <- length(x_coord)-1
    
    
    for (polNb in 1:nbPolygons[[cellNb]]) {
      
      polygon(x=c(x_coord[polNb],x_coord[polNb],x_coord[polNb+1],x_coord[polNb+1]), y=c((countSC-1),countSC,countSC,(countSC-1)), col=colPolygon(x_coord[polNb]), border = NA)  
    }
    countSC <- countSC+1
  }
  dev.off()
  
  #######################
  ##
  ## END OF PLOT CHROM CHUNKS
  ##
  #######################
} else { #end of if there are more than 10 parental cells
  print("No parental cells identified.")
  print("End of the algorithm.")
}  