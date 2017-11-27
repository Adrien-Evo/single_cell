#including functions
#pathToScratch <- "/.mounts/labs/awadallalab/scratch/ialves/scriptsPhasing"
pathToScratch <- getwd()

source(paste(pathToScratch,"/functions_Phasing.R", sep=""))
# .libPaths(c("/.mounts/labs/awadallalab/private/flamaze/R_packages", .libPaths()))
# library(Rmpi)
library(parallel)



#chromosome
chr <- "chr15"
#ind ID
indId <- "AD393"
#step specificity
stepNb <- "phaseTwo"
#computing nb of cells supporting HQ links!! ATENTION: TRUE takes long time
compLinksPhaseOne <- FALSE
exportMatPhaseONe <- FALSE
#computing nb of cells supporting LQ links takes long!!! ATENTION!
compLinksPhaseTwo <- FALSE

##-----------

#FILE NAMES
prefox <- paste(pathToScratch, "/", indId, ".genotypeMatrix.", chr, ".", sep="")
#list of heterozygous SNPs overlapping between genotyping and SC vcf file
hetSNPsList <- paste(pathToScratch, "/", indId, ".hetCallsHQ.", chr, ".c.gt", sep="")
#computation output
output <- paste(pathToScratch, "/",indId, ".", chr, ".phasing.phaseOne.out", sep="")
cat("", file=output, sep="")
#output Matrix for further analysis
linksOutput <- paste(pathToScratch, "/",indId, ".",chr, ".pairwiseComb.phaseOne.txt", sep="")
cat("", file=linksOutput, sep = "")
#output matrix with all the haplotypes
fNameMatrixSortedHap <- paste(pathToScratch, "/",indId, ".", chr, ".hapMatrix_phaseOne.txt", sep="")
#output all pairwise combinations phase two
phaseTwoLinksOutput <- paste(pathToScratch, "/",indId, ".",chr, ".pairwiseComb.phaseTwo.txt", sep="")
#log file for phasing LQ variants
outputLQPhasing <- paste(pathToScratch, "/",indId, ".",chr, ".LQphasing.log", sep="")
cat("", file=outputLQPhasing, sep = "")
##-----------

#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
openGeno <- data.frame(read.table(paste(prefox, "geno", sep=""), header = F, na.strings = "."))
openGenoq <- data.frame(read.table(paste(prefox, "genoq", sep=""), header = F, na.strings = "."))
openHetSNPs <- data.frame(read.table(hetSNPsList, header=T))
namesHapPos <- scan(file=fNameMatrixSortedHap, what = numeric(), nlines = 1)
matrixSortedHap <- data.frame(read.table(fNameMatrixSortedHap, header=T))
names(matrixSortedHap) <- namesHapPos

#keeping positions column
hetSNPsPos <- openHetSNPs[,2]
varNames <- openGeno[,2]
openGeno <- openGeno[,-c(1,2)]
openGenoq <- openGenoq[,-c(1,2)]
# openGeno <- openGeno[,-c(1,2)]
# openGenoq <- openGenoq[,-c(1,2)]

###################
# Quality filter parameters
###################

minHQCells <- 35
minNbCells <- 25
minNbLinks <- 10
##----------

#######################
##
## Compute all pairwise comb for PHASE ONE OR TWO
##
#######################

if (stepNb == "phaseOne") {

  if (compLinksPhaseOne == T) {

    lengthHQ <- apply(openGenoq, 1, function(x) { length(which(x > 20))})
    varNames <- varNames[which(lengthHQ >= minHQCells)] #replace 25 by minHQCells
    HQ_m <- openGeno[which(lengthHQ >= minHQCells),] #replace 25 by minHQCells
    
    #creating a file with the links and corresponding counts
    #this file "linksOutput" is required for step 2 below
    subsettingByQuality_computingHQLinks(HQ_m)

  }
} else if (stepNb == "phaseTwo") {
  
  if (compLinksPhaseTwo == T) {
  
    subsettingByQuality_computingLQLinks(phaseTwoLinksOutput)
    }
} 
####--------------------
##---------------

if (stepNb == "phaseOne" | stepNb == "phaseTwo") {

  if (stepNb == "phaseOne") {
    
      ########################
      ##
      ## STEP 2: import all pairwise comb counts and keep only those with link counts > minNbLinks (GLOBALvariable)
      ##    -----------> CREATE A FUNCTION FOR THIS in order to accommodate into a single script
      ########################
      filePairCom <- "/.mounts/labs/awadallalab/scratch/ialves/scriptsPhasing/AD393.chr15.pairwiseComb.x.txt" #linksOutput
      headerF <- FALSE #false if phase one, true if phase two
      chrNbInFile <- TRUE
      
      PC_tbl <- openLinkCountTable(headerF,chrNbInFile, filePairCom, " ") #change space
      #keep pwise Comb with link counts > minNbLinks
      linksM <- computePCshownByMoreThanXCells(PC_tbl[3:6])
      
      vvv <- which(unlist(lapply(linksM, FUN=length)) == 2) #vector containing those positions with two most prevalent links
      
      #write.table(kkk[vvv,], file="clean_AD393.chr15.pairwiseComb.phase2.txt", quote = F, row.names = F, col.names = F)
      
      #subsetting the original table with pwise Comb counts
      matrix_counts_links_ltTWO <- PC_tbl[vvv,]
      
      #computing a table with the link count order 2,3,4,1 means that the hights nb of counts if in column 2 and 3
      order_m <- t(apply(matrix_counts_links_ltTWO[,3:6],1, function(o) { order(o, decreasing = T) })) #order of the most prevalent links among sites in vvv
      #getting the row indexes 
      colHet <- countingHetLinks(order_m)
      
      ########################
      ##
      ## STEP 3: subset the link count matrix, generate a matrix with all haplotypes and from this a list. This first list contains overlapping haps
      ## change the name of 
      ########################
      #building the m_hap matrix, it includes all the haplotypes that and needs to be clean to rm duplicates and merge overlapping haplotypes
      subMatrix <- matrix_counts_links_ltTWO[colHet,]
      #match(sort(subMatrix[1,3:6], decreasing=T), subMatrix[1,3:6])
      subOrder <- order_m[colHet,]
      
      #creating a matrix with all haplotypes 
      m_hap <- creatingHapMatrixFromLinkCountMatrix(subMatrix,subOrder)
      
      #transforming m_hap into a list
      l_hap <- creatingHapListFromLinkCountMatrix(m_hap)
      
      
      #the next steps remove direct duplicates followed by merging them (and merging partial duplicates - sharing the same var Names but with the inverted
      #haplotypes)
      uniqueHap <- unique(l_hap)
      #generating a list of unique merged haplotypes: FINAL config
      tmpListToRm <- excludingDuplicates(uniqueHap)
      hapList <- tmpListToRm[[1]]
      namesHapList <- tmpListToRm[[2]]
      rm(tmpListToRm)
      
      if (exportMatPhaseONe == T) {
        
        #converting hapList into matrix to EXPORT
        finalHapMatrix <- convertingIntoMatrix(hapList,namesHapList)
        
        #save final haplotype matrix
        write.table(finalHapMatrix, file=fNameMatrixSortedHap, quote = F, row.names = F, col.names = T, sep="\t")
        
      }
      
  } else if (stepNb == "phaseTwo") {
      
    ########################
    ##
    ## STEP 2: import all pairwise comb counts and keep only those with link counts > minNbLinks (GLOBALvariable)
    ##    -----------> CREATE A FUNCTION FOR THIS in order to accommodate into a single script
    ########################
    
    filePairCom <- paste(pathToScratch,"/AD393.chr15.pairwiseComb.x.txt", sep="") #linksOutput
    headerF <- FALSE #false if phase one, true if phase two
    chrNbInFile <- TRUE
    
    PC_tbl <- openLinkCountTable(headerF,chrNbInFile, filePairCom, " ") #change space
    #keep pwise Comb with link counts > minNbLinks
    linksM <- computePCshownByMoreThanXCells(PC_tbl[3:6])
    
    vvv <- which(unlist(lapply(linksM, FUN=length)) == 2) #vector containing those positions with two most prevalent links
    
    #write.table(kkk[vvv,], file="clean_AD393.chr15.pairwiseComb.phase2.txt", quote = F, row.names = F, col.names = F)
    
    #subsetting the original table with pwise Comb counts
    matrix_counts_links_ltTWO <- PC_tbl[vvv,]
    
    #computing a table with the link count order 2,3,4,1 means that the hights nb of counts if in column 2 and 3
    order_m <- t(apply(matrix_counts_links_ltTWO[,3:6],1, function(o) { order(o, decreasing = T) })) #order of the most prevalent links among sites in vvv
    #getting the row indexes 
    colHet <- countingHetLinks(order_m)
    
    ########################
    ##
    ## STEP 3: subset the link count matrix, generate a matrix with all haplotypes and from this a list. This first list contains overlapping haps
    ## change the name of 
    ########################
    #building the m_hap matrix, it includes all the haplotypes that and needs to be clean to rm duplicates and merge overlapping haplotypes
    subMatrix <- matrix_counts_links_ltTWO[colHet,]
    #match(sort(subMatrix[1,3:6], decreasing=T), subMatrix[1,3:6])
    subOrder <- order_m[colHet,]
    
    #creating a matrix with all haplotypes 
    m_hap <- creatingHapMatrixFromLinkCountMatrix(subMatrix,subOrder)
    
    #transforming m_hap into a list
    l_hap <- creatingHapListFromLinkCountMatrix(m_hap)
    
    
    #the next steps remove direct duplicates followed by merging them (and merging partial duplicates - sharing the same var Names but with the inverted
    #haplotypes)
    uniqueHap <- unique(l_hap)
    #generating a list of unique merged haplotypes: FINAL config
    tmpListToRm <- excludingDuplicates(uniqueHap)
    hapList <- tmpListToRm[[1]]
    rm(tmpListToRm)
    rm(uniqueHap)
    #####--------------------
    ##-------------
    #-------
    
    #resetting QC filters
    minHQCells <- 25
    minNbCells <- 4
    minNbLinks <- 2 
    
    ########################
    ##
    ## STEP 2: import all pairwise comb counts and keep only those with link counts > minNbLinks (GLOBALvariable)
    ##  CREATE A FUNCTION FOR THIS in order to accommodate into a single scripts
    ########################
    filePairCom <- paste(pathToScratch, "/AD393.chr15.pairwiseComb.phase2.x.txt", sep="")
    headerF <- TRUE #false if phase one, true if phase two
    chrNbInFile <- FALSE
    
    PC_tbl <- openLinkCountTable(headerF,chrNbInFile, filePairCom, "\t") #change space
    #keep pwise Comb with link counts > minNbLinks
    linksM <- computePCshownByMoreThanXCells(PC_tbl[3:6])
    
    vvv <- which(unlist(lapply(linksM, FUN=length)) == 2) #vector containing those positions with two most prevalent links
    
    #write.table(kkk[vvv,], file="clean_AD393.chr15.pairwiseComb.phase2.txt", quote = F, row.names = F, col.names = F)
    
    #subsetting the original table with pwise Comb counts
    matrix_counts_links_ltTWO <- PC_tbl[vvv,]
    
    #computing a table with the link count order 2,3,4,1 means that the hights nb of counts if in column 2 and 3
    order_m <- t(apply(matrix_counts_links_ltTWO[,3:6],1, function(o) { order(o, decreasing = T) })) #order of the most prevalent links among sites in vvv
    #getting the row indexes 
    colHet <- countingHetLinks(order_m)
    
    ########################
    ##
    ## STEP 3: subset the link count matrix, generate a matrix with all haplotypes and from this a list. This first list contains overlapping haps
    ##
    ########################
    subMatrix <- matrix_counts_links_ltTWO[colHet,]
    subOrder <- order_m[colHet,]
    nrow(subMatrix)
    names_subM <- sort(union(subMatrix[,1],subMatrix[,2]))
    varLQstep <- names(table(subMatrix[,1]))
    finalHapList <- list()
    reducedFinalHapList <- list()
    
    # library(parallel)
    # 
    # # Calculate the number of cores
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    
    # Initiate cluster
    #cl <- makeCluster( mpi.universe.size(), type = 'MPI')
    clusterExport(cl, c("hapList", "subMatrix", "subOrder","outputLQPhasing","creatingLQhaplotypes", "excludingDuplicates", "mergingHaplotypes", "namesHapListFunction"))
    finalHapList <- parLapply(cl=cl, varLQstep[1:length(varLQstep)], creatingLQhaplotypes)
    #finalHapList
    reducedFinalHapList <- list()

    countLQhap <- 1

    for (pos in 1:length(finalHapList)) {

     if (!is.null(colnames(finalHapList[[pos]][[1]]))) {

       reducedFinalHapList[[countLQhap]] <- finalHapList[[pos]][[1]]
       countLQhap <- countLQhap+1
      }

   }
   rm(finalHapList)
   mergedHapList <- mergingHaplotypes(reducedFinalHapList)
   rm(reducedFinalHapList)
  }

}



namesPhased <- unlist(lapply(mergedHapList, function(x) {colnames(x)}))
indexHap <- c()
for (hap_index in 1:length(mergedHapList)) {
 indexHap <- c(indexHap, rep(hap_index, length(colnames(mergedHapList[[hap_index]]))))

}
resHapPhased <- excludingDuplicates(mergedHapList)
mHap <- matrix(unlist(resHapPhased[[1]]), ncol=length(colnames(resHapPhased[[1]][[1]])))
colnames(mHap) <- colnames(resHapPhased[[1]][[1]])
write.table(mHap, file=paste(pathToScratch, "/1158Var_haplotypes.txt", sep=""), quote = F, row.names = F, col.names = T, sep="\t")


stopCluster(cl)
#mpi.exit()

