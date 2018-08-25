#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Sourcing funtions #TODO package it

source("./computingLinks_functions.v2.R")
source("./subsettingLinksMatrix_function.v2.R")
source("./list_to_matrix_to_list_function.v2.R")
source("./creatingHapList_gettingHapNames_functions.v2.R")
source("./excludingDuplicates_function.v2.R")
source("./mergingHap_function.v2.R")
source("./creatingHapByMerging.v3.R")
#library(Rmpi)
library(parallel)
dateToday <- Sys.Date()
nbOfCores <- 10L

#chromosome
chr <- args[2]
#ind ID
indId <- args[1]

##Parameters
#HQ ratio difference to accept link
HQ_ratio <- 5
#LQ ratio difference to accept link
LQ_ratio <- 3
#proportion of closer HQhaps
propHQhaps <- 1
#export HQ haplotypes
exportMatPhaseONe <- T
#export HQ links table
exportHQLinksTbl <- F
#export LQ links table
exportLQLinksTbl <- F
#export HQ link table
printHQcountM <- T
#nb of LQ links supporting a HQ_HQ connection
minLQsupport <- 1 #this is just in case of using creatingHQ_LQvar links with tags
##-----------
###################
# Quality filter parameters
###################

minHQCells <- 25
minNbCells <- 10
minNbLinks <- 5
##----------

##Names of outputs files
folderName <- paste(pathToScratch, "/",indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, "_", dateToday, sep = "")

##This will be removed since you have function, to write temp files
dir.create(folderName)
pathToTmp <- paste(folderName,"tmp", sep="/") #added by March 2
dir.create(pathToTmp) #added by March 2

#FILE NAMES
prefox <- paste0("/.mounts/labs/awadallalab/scratch/SingleCell/vcfs_sc/allSamples", "/", indId, "/", indId, ".SingleCellsHetSNPs.", chr)
#log file containing info on the main function
mainOutput <- paste(folderName, "/", indId, ".", chr, ".mainLogFile_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".log", sep="")
#Cat is there to create the file *this is bad practise
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
#final LQ variants anchoring HQ haplotypes
LQ_HQ_matrixFile <- paste(folderName, "/", indId, ".",chr, ".LQvar_HQhap_matrix_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks,".txt", sep="")

##-----------

#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
##All the cat in main output are log files that should be replace by futile logging R pckg
cat("Opening genomic matrices...", file=mainOutput, sep="\n", append = T)
openGeno <- data.frame(read.table(paste0(prefox, ".GT.FORMAT"), header = T, na.strings = "."))
openGenoq <- data.frame(read.table(paste0(prefox, ".GQ.FORMAT"), header = T, na.strings = "."))

#keeping positions column
#hetSNPsPos <- openHetSNPs[,2]
#Keep the positions in a separate vector
varNames <- openGeno[,2]
#openGeno <- openGeno[,-c(1,2)]
#openGenoq <- openGenoq[,-c(1,2)]
#Settings some parameters from the files
nbOfCellCovPerSite <- unlist(lapply(1:nrow(openGeno), function(x) { sum(!is.na(openGeno[x,])) }))
#hist(nbOfCellCovPerSite)
nbOfSitesPerCell <- unlist(lapply(3:ncol(openGeno), function(x) { sum(!is.na(openGeno[,x])) }))
# #hist(nbOfSitesPerCell)
quant95 <- quantile(nbOfSitesPerCell, probs=c(0.025, 0.85)) #changed by March15
# #hist(nbOfSitesPerCell[which(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2])])
#Removing outliers cells lmooking at the two above parameters
#Threshold are hardcoded : needs to be change in parameters
rmCells <- which(!(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2]))
openGeno <- openGeno[,-rmCells]
openGenoq <- openGenoq[,-rmCells]

#Removing the chr and postion
openGeno <- openGeno[,-c(1,2)]
openGenoq <- openGenoq[,-c(1,2)]

#garbage collect but not sure if its useful here yet
gc()
#######################
##
## Compute all pairwise comb for PHASE ONE 
##
#######################
#LOGGING not in futile hqs to be changed
cat("Retrieving HQ variants.", file=mainOutput, sep="\n", append = T)
cat("HQvar-HQvar will be done with the following QC: ", file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells with variant genotyped at GQ > 20: ", minHQCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each HQvar-HQvar combination: ", minNbCells, ".", sep = ""), file=mainOutput, sep="\n", append = T)
cat(paste("Number of cells covering each pairwise combination of alleles: ", minNbLinks, ".", sep = ""), file=mainOutput, sep="\n", append = T)

##Checking for each positions the number of cells with a genotype quality bigger than 20
lengthHQ <- apply(openGenoq, 1, function(x) { length(which(x > 20))})
##Storing the sites with more than minHQcells across all cells
HQvarNames <- varNames[which(lengthHQ >= minHQCells)] #replace 25 by minHQCells

#COMMENT
cat(paste("Total number of HQ variants: ", length(HQvarNames), ".", sep = ""), file=mainOutput, sep="\n", append = T)

#Subset of GT matrix satisfying the criteria of cells with good HQ
#Maybe the quality matrix is not needed at this anymore
HQ_m <- openGeno[which(lengthHQ >= minHQCells),] #replace 25 by minHQCells

#COMMENT
cat("Creating cluster for HQ phasing...", file=mainOutput, sep="\n", append = T)
links_and_order_tmp <- list()
#creating a file with the links and corresponding counts
#this file "linksOutput" is required for step 2 below
#computingLinks <- function(tmpSite, listOfVar, HQ_genotypeMatrix, outFileName)
##First occurence of the function computing links
##HQVARNAMES is actuall yt HQSITES shoudl be renamed !!!!! 
links_and_order_tmp <- mclapply(HQvarNames, function(x) { computingLinks(tmpSite=x, listOfVar=HQvarNames, HQ_genotypeMatrix=HQ_m, link_r = HQ_ratio, outFileName=outputPhaseOne) }, mc.cores = nbOfCores) #changed April 19
#COMMENT
cat("HQ phasing... DONE.", file=mainOutput, sep="\n", append = T)

# listLinks <- listLinks[which(sapply(listLinks, length) > 1)] #modified by Nov 27
# 
# 
# links_and_order_tmp <- mclapply(listLinks, function(x) {  subsettingLinksMatrix(x, HQ_ratio)}, mc.cores = nbOfCores)

# fullHQlinksComb <- data.frame()
# fullHQlinksComb <- do.call(rbind, listLinks[which(sapply(listLinks, length) > 1)]) #removed by Nov 27

links_and_order_tmp <- links_and_order_tmp[!sapply(links_and_order_tmp, is.null)] #modified by Nov 27
#---

##First use of one the auxilary function
l_hap <- mclapply(links_and_order_tmp, FUN = creatingHQhapList, mc.cores = nbOfCores)

if (length(l_hap) > 1) { #added by March 7
  #generating a list of unique merged haplotypes: FINAL config
  tmpListToRm <- excludingDuplicates(l_hap)
} else { #added by March 7
  tmpListToRm <- l_hap
}

##Check when and how this stuff is used later
hapList <- tmpListToRm
#COMMENT
cat(paste("HQ haplotypes list contains:", length(hapList), "HQ haplotypes.", sep = " "), file=mainOutput, sep="\n", append = T)

##Function auxilary
namesHapList <- namesHapListFunction(hapList)
rm(links_and_order_tmp)
rm(tmpListToRm)
rm(l_hap)
##Getting the trash out the way
gc()

##It's an option in the beginning
if (exportMatPhaseONe == T) {
  #COMMENT
  cat("Exporting matrix with all HQ haplotypes.", file=mainOutput, sep="\n", append = T)
  #converting hapList into matrix to EXPORT
  ##AUx function text formatting. Not algo related to bio
  finalHapMatrix <- convertingIntoMatrix(hapList,namesHapList)
  finalHapMatrix <- finalHapMatrix[,order(as.numeric(colnames(finalHapMatrix)))]
  
  ##Need to check that. The ouput of the function is not always the same so thats why the matrix format neeeds to
  ##be forced. Can be a problem TODO
  
  if (is.matrix(finalHapMatrix)) {
    
    #save final haplotype matrix
    write.table(finalHapMatrix, file=fNameMatrixSortedHap, quote = F, row.names = F, col.names = T, sep="\t")
    
  } else {
    ##
    write.table(t(as.matrix(finalHapMatrix)), file=fNameMatrixSortedHap, quote = F, row.names = F, col.names = T, sep="\t")
    
    
  }
  
}

#######################
##
## Compute all pairwise comb between LQ var and HQ var
##Second par tof the process. takes the longest time. Maybe test that afterwards
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
varOutHQhaps <- setdiff(varNames, sort(as.numeric(varInHQhaps), decreasing = F)) #added by Apr 4
rm(namesHapList)
#COMMENT
cat("Creating cluster for LQ-HQhap links computation...", file=mainOutput, sep="\n", append = T)

#initializing cluster
#LQ_HQ_links_list <- list()
links_and_order_List <- list()
#openGenoq, openGeno, varNames, outputLQPhasing, minHQCells
#changed by Apr 4
links_and_order_List <- mclapply(varOutHQhaps, function(x) { subsettingByQuality_computingLQLinks(tmpLQSite=x, varNamesInHQhaps=varInHQhaps, ratio=LQ_ratio, outLog=LQlinksLog) }, mc.cores = nbOfCores)
#change from varNames to varNames-HQvar
#COMMENT
cat("LQvar-HQhap links computation... DONE.", file=mainOutput, sep="\n", append = T)

#-- NEW OCT 25
#-- CHANGED APRIL 4
#COMMENT
#cat("Creating cluster for subsetting LQ-HQhap links per LQvar matrices...", file=mainOutput, sep="\n", append = T)
# links_and_order_List <- list()
# links_and_order_List <- mclapply(LQ_HQ_links_list, function(x) {  subsettingLinksMatrix(x, LQ_ratio)}, mc.cores = nbOfCores)
#COMMENT
cat("Subsetting and ordering of LQ-HQhap links per LQvar matrices... DONE.", file=mainOutput, sep="\n", append = T)
links_and_order_List <- links_and_order_List[!sapply(links_and_order_List, is.null)] #modified by Oct 26
cat(paste0("Links and order list size: ", length(links_and_order_List), "."), file=mainOutput, sep="\n", append = T) #added by late Nov 27
cat(paste0("Links and order list size in bytes is: ", object.size(links_and_order_List), "."), file=mainOutput, sep="\n", append = T) #added by Jan 24

#---
if (exportLQLinksTbl) {
  
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
  
}
rm(LQ_HQ_links_list) #modified by Oct 25 
gc()
# ##########################
# ##
# ##Step building the final matrix with recomb map. HQ + LQ . A bit like steps 2b
# ###########################
#COMMENT
cat("Building up LQvar-HQhap connection.", file=mainOutput, sep="\n", append = T)

#computing the abs number of closest haplotypes           for LQ regarding HQ. Want to look at closer haplotypoes to decrease the chance of recomb break
nbOfCloserHQhap <- round(length(hapList)*propHQhaps, digits = 0)
if (nbOfCloserHQhap == 0) { #in case the nb of HQhaps is too small and/or the chosen prop.
  nbOfCloserHQhap <- 1
}

#COMMENT
cat("Creating cluster for LQvar-HQhap phasing...", file=mainOutput, sep="\n", append = T)
## Calculate the number of cores We dont kknow what it means
cat(paste("LQVar", paste0("HQ_", 1:length(hapList), collapse = "\t"), sep = "\t"), file=LQ_HQ_matrixFile, sep="\n")

s_list <- seq(1, length(links_and_order_List), 10)
for (ind_l in s_list) {
  
  tmpListToRm <- mclapply(links_and_order_List[c(ind_l:(ind_l+9))], function(varPos) { creatingLQhaplotypes(varPos) }, mc.cores = nbOfCores)

}  
  cat("Computation of phasing done! Printing LQvar HQhap matrix!", file=mainOutput, sep="\n", append = T)

setwd(pathToTmp)
system(paste0("cat *.txt > ", LQ_HQ_matrixFile))
setwd(folderName)
unlink(pathToTmp, recursive = T)
cat("LQvar HQhap matrix written. Tmp files deleted. All good.", file=mainOutput, sep="\n", append = T)


#####################
##
## End of the main function. Beautiful. Infinite Bisous
##
#####################
