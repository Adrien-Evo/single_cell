#!/usr/bin/env Rscript

#including functions
pathToScratch <- "/.mounts/labs/awadallalab/scratch/ialves/scriptsPhasing"
#pathToScratch <- "/Users/isabelalves/Documents/OICR/singleCells/"

source(paste(pathToScratch,"/functions_new_Oct27/computingLinks_functions.R", sep=""))
source(paste(pathToScratch,"/functions_new_Oct27/subsettingLinksMatrix_function.R", sep=""))
<<<<<<< HEAD
source(paste(pathToScratch,"/functions_new_Oct27/list_to_matrix_to_list_function.v1.R", sep=""))
=======
source(paste(pathToScratch,"/functions_new_Oct27/list_to_matrix_to_list_function.R", sep=""))
>>>>>>> 7f3eaabb039066a8eac7b49e27c5f0ac6afdbf16
source(paste(pathToScratch,"/functions_new_Oct27/creatingHapList_gettingHapNames_functions.v1.R", sep=""))
source(paste(pathToScratch,"/functions_new_Oct27/excludingDuplicates_function.v1.R", sep=""))
source(paste(pathToScratch,"/functions_new_Oct27/mergingHap_function.R", sep=""))
source(paste(pathToScratch,"/functions_new_Oct27/creatingHapByMerging.v1.R", sep=""))
.libPaths( c("/.mounts/labs/awadallalab/private/flamaze/R_packages", .libPaths() ) )
#library(Rmpi)
library(parallel)
dateToday <- "TODAY"
nbOfCores <- 10L

#chromosome
chr <- "CHROM"
#ind ID
indId <- "INDIVIDUAL"
#HQ ratio difference to accept link
HQ_ratio <- 5
#LQ ratio difference to accept link
LQ_ratio <- 3
#proportion of closer HQhaps
propHQhaps <- 0.10
#export HQ haplotypes
exportMatPhaseONe <- T
#export HQ links table
exportHQLinksTbl <- F
#export LQ links table
exportLQLinksTbl <- F
#export HQ link table
printHQcountM <- T
#rm outlier cells
rmBadCells <- T
#nb of LQ links supporting a HQ_HQ connection
minLQsupport <- 1 #this is just in case of using creatingHQ_LQvar links with tags
##-----------
###################
# Quality filter parameters
###################

minHQCells <- 25
minNbCells <- 20
minNbLinks <- 10
##----------

folderName <- paste(pathToScratch, "/",indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, "_", dateToday, sep = "")
dir.create(folderName)
#FILE NAMES
prefox <- paste(pathToScratch, "/", indId, ".genotypeMatrix.", chr, ".", sep="")
#list of heterozygous SNPs overlapping between genotyping and SC vcf file
hetSNPsList <- paste(pathToScratch, "/", indId, ".hetCallsHQ.", chr, ".c.gt", sep="")
#log file containing info on the main function
mainOutput <- paste(folderName, "/", indId, ".", chr, ".mainLogFile_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".log", sep="")
cat("", file=mainOutput, sep="")
#HQvar_HQvar into HQhaps computation output
outputPhaseOne <- paste(folderName, "/", indId, ".", chr, ".HQ_HQvarLinks", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".log", sep="")
cat("", file=outputPhaseOne, sep="")
#output Matrix for further analysis
linksOutput <- paste(folderName, "/", indId, ".",chr, ".phaseOneHQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".txt", sep="")
#output matrix with all the haplotypes
fNameMatrixSortedHap <- paste(folderName, "/", indId, ".", chr, ".HQhapMatrix", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".txt", sep="")
#output all pairwise combinations phase two
phaseTwoLinksOutput <- paste(folderName, "/", indId, ".",chr, ".phaseTwoHQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks,".txt", sep="")
#log file for phasing LQ variants
outputLQPhasing <- paste(folderName, "/", indId, ".",chr, ".finalPhasing_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".log", sep="")
cat("", file=outputLQPhasing, sep = "")
#log file for link computation with LQ
LQlinksLog <- paste(folderName, "/", indId, ".",chr, ".LQvar_HQhapLinks", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".log", sep="")
cat("", file=LQlinksLog, sep = "")
#final haplotype file
finalCompleteHapFile <- paste(folderName, "/", indId, ".",chr, ".finalHapNeighborHaps_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks,".txt", sep="")
##-----------

#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
cat("Opening genomic matrices...", file=mainOutput, sep="\n", append = T)
openGeno <- data.frame(read.table(paste(prefox, "geno", sep=""), header = F, na.strings = "."))
openGenoq <- data.frame(read.table(paste(prefox, "genoq", sep=""), header = F, na.strings = "."))
#openHetSNPs <- data.frame(read.table(hetSNPsList, header=T))
# namesHapPos <- scan(file=fNameMatrixSortedHap, what = numeric(), nlines = 1)
# matrixSortedHap <- data.frame(read.table(fNameMatrixSortedHap, header=T))
# names(matrixSortedHap) <- namesHapPos

#keeping positions column
varNames <- openGeno[,2]

if (rmBadCells) { #added by Dec 2
  
  nbOfCellCovPerSite <- unlist(lapply(1:nrow(openGeno), function(x) { sum(!is.na(openGeno[x,])) }))
  #hist(nbOfCellCovPerSite)
  nbOfSitesPerCell <- unlist(lapply(3:ncol(openGeno), function(x) { sum(!is.na(openGeno[,x])) }))
  # #hist(nbOfSitesPerCell)
  quant95 <- quantile(nbOfSitesPerCell, probs=c(0.025, 0.975))
  # #hist(nbOfSitesPerCell[which(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2])])
  rmCells <- which(!(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2]))
  openGeno <- openGeno[,-rmCells]
  openGenoq <- openGenoq[,-rmCells]
}  #added by Dec 2
openGeno <- openGeno[,-c(1,2)]
openGenoq <- openGenoq[,-c(1,2)]

#######################
##
## Compute all pairwise comb for PHASE ONE 
##
#######################
#COMMENT
cat("Retrieving HQ variants.", file=mainOutput, sep="\n", append = T)
cat("HQvar-HQvar will be done with the following QC: ", file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells with variant genotyped at GQ > 20: ", minHQCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each HQvar-HQvar combination: ", minNbCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each pairwise combination of alleles: ", minNbLinks, ".", sep = ""), file=mainOutput, sep="\n", append = T)

lengthHQ <- apply(openGenoq, 1, function(x) { length(which(x > 20))})
HQvarNames <- varNames[which(lengthHQ >= minHQCells)] #replace 25 by minHQCells

#COMMENT
cat(paste("Total number of HQ variants: ", length(HQvarNames), ".", sep = ""), file=mainOutput, sep="\n", append = T)

HQ_m <- openGeno[which(lengthHQ >= minHQCells),] #replace 25 by minHQCells

#COMMENT
cat("Creating cluster for HQ phasing...", file=mainOutput, sep="\n", append = T)

listLinks <- list()
#creating a file with the links and corresponding counts
#this file "linksOutput" is required for step 2 below
#computingLinks <- function(tmpSite, listOfVar, HQ_genotypeMatrix, outFileName) 
listLinks <- mclapply(HQvarNames, function(x) { computingLinks(tmpSite=x, listOfVar=HQvarNames, HQ_genotypeMatrix=HQ_m, outFileName=outputPhaseOne) }, mc.cores = nbOfCores)
#COMMENT
cat("HQ phasing... DONE.", file=mainOutput, sep="\n", append = T)

listLinks <- listLinks[which(sapply(listLinks, length) > 1)] #modified by Nov 27

links_and_order_tmp <- list()
links_and_order_tmp <- mclapply(listLinks, function(x) {  subsettingLinksMatrix(x, HQ_ratio)}, mc.cores = nbOfCores)

# fullHQlinksComb <- data.frame()
# fullHQlinksComb <- do.call(rbind, listLinks[which(sapply(listLinks, length) > 1)]) #removed by Nov 27

links_and_order_tmp <- links_and_order_tmp[!sapply(links_and_order_tmp, is.null)] #modified by Nov 27
#---

#merge matrices over all the LQ var into a single matrix
# subMatrix <- do.call(rbind, (lapply(links_and_order_tmp, "[[", 1))) #modified by Nov 27
# subOrder <- do.call(rbind, (lapply(links_and_order_tmp, "[[", 2))) #modified by Nov 27 



# #links_and_order_tmp <- subsettingLinksMatrix(fullHQlinksComb, HQ_ratio) #removed by Nov 27
# subMatrix <- links_and_order_tmp[[1]] #removed by Nov 27
# subOrder <- links_and_order_tmp[[2]] #removed by Nov 27

# if (printHQcountM) {
#   #COMMENT
#   cat("Exporting table with HQ-HQvar links.", file=mainOutput, sep="\n", append = T)
#   write.table(subMatrix, file=linksOutput, quote = F, row.names = F, col.names = T, sep="\t", append = T)
#   
# }



#rm(fullHQlinksComb) #removed by Nov 27
#creating a matrix with all haplotypes 
#m_hap <- creatingHapMatrixFromLinkCountMatrix(subMatrix,subOrder) #uncommented by Dec 1

#transforming m_hap into a list
#l_hap <- creatingHapListFromLinkCountMatrix(m_hap)  #uncommented by Dec 1

l_hap <- mclapply(links_and_order_tmp, FUN = creatingHQhapList, mc.cores = nbOfCores)

#### COMMENT OUT if you are generating the hapList
# scanNames <- scan(fNameMatrixSortedHap, nlines = 1,what = character(), sep = "\t")
# tmpM <- read.table(fNameMatrixSortedHap, header=T)
# colnames(tmpM) <- scanNames
#l_hap <- creatingHapListFromLinkCountMatrix(tmpM)
###----------
#---

#the next steps remove direct duplicates followed by merging them (and merging partial duplicates - sharing the same var Names but with the inverted
#haplotypes)
uniqueHap <-  mclapply(l_hap, FUN = creatingHapList, mc.cores = nbOfCores)  # changed by Dec 1
uniqueHap <- unique(uniqueHap)
#generating a list of unique merged haplotypes: FINAL config
tmpListToRm <- excludingDuplicates(uniqueHap)
tmpListToRm <-  mclapply(tmpListToRm, FUN = creatingHapList, mc.cores = nbOfCores)  # changed by Dec 1
#rm(uniqueHap)


hapList <- tmpListToRm
#COMMENT
cat(paste("HQ haplotypes list contains:", length(hapList), "HQ haplotypes.", sep = " "), file=mainOutput, sep="\n", append = T)

namesHapList <- namesHapListFunction(hapList)
rm(links_and_order_tmp)
rm(tmpListToRm)
# rm(subMatrix)
# rm(subOrder)
rm(uniqueHap)
#rm(m_hap)
rm(l_hap)

#convert all elements of the hapList into the same type of objects
<<<<<<< HEAD
hapList <- mclapply(hapList, FUN = creatingHapList, mc.cores = nbOfCores) # changed by Dec 1 PROBABLY REMOVE!!!!
=======
hapList <- mclapply(hapList, FUN = creatingHapList, mc.cores = nbOfCores) # changed by Dec 1
>>>>>>> 7f3eaabb039066a8eac7b49e27c5f0ac6afdbf16

if (exportMatPhaseONe == T) {
  #COMMENT
  cat("Exporting matrix with all HQ haplotypes.", file=mainOutput, sep="\n", append = T)
  #converting hapList into matrix to EXPORT
  finalHapMatrix <- convertingIntoMatrix(hapList,namesHapList)
  finalHapMatrix <- finalHapMatrix[,order(as.numeric(colnames(finalHapMatrix)))]
  #save final haplotype matrix
  write.table(finalHapMatrix, file=fNameMatrixSortedHap, quote = F, row.names = F, col.names = T, sep="\t")
  
}

#######################
##
## Compute all pairwise comb between LQ var and HQ var
##
#######################
#COMMENT
# cat("LQ variants are required to link HQ haplotypes.", file=mainOutput, sep="\n", append = T)
#resetting QC filters
minHQCells <- 10
minNbCells <- 4
minNbLinks <- 2 

#COMMENT
cat("LQvar-HQhaplotype links will be done with the following QC: ", file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells with variant genotyped at GQ > 20: ", minHQCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each LQvar-varWithinHQhap combination: ", minNbCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each LQvar-varWithinHQhap pairwise combination of alleles: ", minNbLinks, ".", sep = ""), file=mainOutput, sep="\n", append = T)
# 
varInHQhaps <- namesHapList
rm(namesHapList)
#COMMENT
cat("Creating cluster for LQ-HQhap links compution...", file=mainOutput, sep="\n", append = T)

#initializing cluster
LQ_HQ_links_list <- list()
#openGenoq, openGeno, varNames, outputLQPhasing, minHQCells
LQ_HQ_links_list <- mclapply(varNames, function(x) { subsettingByQuality_computingLQLinks(tmpLQSite=x, varNamesInHQhaps=varInHQhaps, outLog=LQlinksLog) }, mc.cores = nbOfCores)
#change from varNames to varNames-HQvar
#COMMENT
cat("LQvar-HQhap links computation... DONE.", file=mainOutput, sep="\n", append = T)

#-- NEW OCT 25
#COMMENT
cat("Creating cluster for subsetting LQ-HQhap links per LQvar matrices...", file=mainOutput, sep="\n", append = T)

links_and_order_List <- list()
links_and_order_List <- mclapply(LQ_HQ_links_list, function(x) {  subsettingLinksMatrix(x, LQ_ratio)}, mc.cores = nbOfCores)
#COMMENT
cat("Subsetting and ordering of LQ-HQhap links per LQvar matrices... DONE.", file=mainOutput, sep="\n", append = T)
links_and_order_List <- links_and_order_List[!sapply(links_and_order_List, is.null)] #modified by Oct 26
cat(paste0("Links and order list size: ", length(links_and_order_List), "."), file=mainOutput, sep="\n", append = T) #added by late Nov 27
#---
if (exportLQLinksTbl) {
<<<<<<< HEAD
  
=======

>>>>>>> 7f3eaabb039066a8eac7b49e27c5f0ac6afdbf16
  #merge matrices over all the LQ var into a single matrix
  subMatrix <- do.call(rbind, (lapply(links_and_order_List, "[[", 1))) #modified by Oct 26
  #subOrder <- do.call(rbind, (lapply(links_and_order_List, "[[", 2))) #modified by Oct 26
  
  #### COMMENT OUT if you are generating the phaseTwoLinksOutput matrix
  #fullHQlinksComb <- read.table(phaseTwoLinksOutput, header=T)
  ###--------
  #-----
  
  #COMMENT
  cat("Exporting matrix with supported LQvar-HQvar links.", file=mainOutput, sep="\n", append = T)
  write.table(subMatrix, file=phaseTwoLinksOutput, quote = F, row.names = F, col.names = T, sep="\t")
<<<<<<< HEAD
  
=======

>>>>>>> 7f3eaabb039066a8eac7b49e27c5f0ac6afdbf16
}
rm(LQ_HQ_links_list) #modified by Oct 25 
# ##########################
# ##
# ##Step THREE
# ###########################
#COMMENT
cat("Building up LQvar-HQhap connection.", file=mainOutput, sep="\n", append = T)

#computing the abs number of closer haplotypes
nbOfCloserHQhap <- round(length(hapList)*propHQhaps, digits = 0)
#COMMENT
cat("Creating cluster for LQvar-HQhap phasing...", file=mainOutput, sep="\n", append = T)
## Calculate the number of cores
finalHapList <- list()
finalHapList <- mclapply(links_and_order_List, function(varPos) { creatingLQhaplotypes(varPos) }, mc.cores = nbOfCores)

#COMMENT
cat("LQvar-HQhap phasing... DONE.", file=mainOutput, sep="\n", append = T)
#-
#removing empty entries from the finalHapList
reducedFinalHapList <- list()
reducedFinalHapList <- finalHapList[!sapply(finalHapList, is.null)] #modified by Oct 26
rm(finalHapList)
#organizing the list
## Calculate the number of cores
# no_cores <- detectCores() #modified by Oct 27
# cl <- makeCluster(no_cores) #modified by Oct 27
# clusterExport(cl, c("reducedFinalHapList", "creatingHapList")) #modified by Oct 27
cleanReducedList <- list() #modified by Oct 27
<<<<<<< HEAD
cleanReducedList <- mclapply(reducedFinalHapList, FUN = creatingHapList, mc.cores = nbOfCores) #modified by Dec 1

=======
cleanReducedList <- mclapply(reducedFinalHapList, FUN =  creatingHapList, mc.cores = nbOfCores) #modified by Dec 1
 
>>>>>>> 7f3eaabb039066a8eac7b49e27c5f0ac6afdbf16
#merging the HQhaps with the LQ_HQhap links
mergedHapList <- excludingDuplicates(cleanReducedList)
#COMMENT
cat("Merging all LQvar-HQhaplotypes ... DONE.", file=mainOutput, sep="\n", append = T)

if (length(mergedHapList) > 1) {
  #COMMENT
  cat("Multiple haplotypes still exist. ERROR.", file=mainOutput, sep="\n", append = T)
  
} else if (length(mergedHapList) == 1) {
  #COMMENT
  cat("Parent phasing SUCCESSFUL. Well done!", file=mainOutput, sep="\n", append = T)
  #ordering the final hap
  matrixFinalHap <- mergedHapList[[1]][,match(sort(as.numeric(names(mergedHapList[[1]]))), names(mergedHapList[[1]]))] 
  write.table(matrixFinalHap, file=finalCompleteHapFile, quote = F, row.names = F, col.names = T, sep = "\t")
  
}

#####################
##
## End of the main function
##
#####################
