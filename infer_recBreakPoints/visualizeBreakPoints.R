#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

#pathToScratch <- getwd()
#pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/AD393/scriptsPhasing"
#including functions
pathToScratch <- "/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"

###################
#declaring variables
###################
dateToday <- args[1]
#chromosome
chr <- args[2]
#ind ID
indId <- args[3]
badCellsRmved <- T

###################
# Quality filter parameters
###################

minHQCells <- 25
minNbCells <- 10
minNbLinks <- 5
##----------

folderName <- paste(pathToScratch, "/",indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, "_", dateToday, sep = "")
fileName <- paste(indId, ".", chr, ".finalHapNeighborHaps_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, ".txt", sep = "")
outfileName <- paste(folderName, "/SC_phasedHapHMap_", indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, ".pdf", sep = "")

#opening haplotype file
hapFile <- scan(paste(folderName, "/", fileName, sep = ""), nlines=1, what=numeric())
openHap <- read.table(paste(folderName, "/", fileName, sep = ""), header=T, sep="\t")
colnames(openHap) <- hapFile


#GENO MATRIX FILE NAMES
prefox <- paste0("/.mounts/labs/awadallalab/scratch/SingleCell/vcfs_sc/allSamples", "/", indId, "/", indId, ".SingleCellsHetSNPs.", chr)
#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
openGeno <- data.frame(read.table(paste0(prefox, ".GT.FORMAT"), header = F, na.strings = "."))
#openGenoq <- data.frame(read.table(paste(prefox, "genoq", sep=""), header = F, na.strings = "."))

if (badCellsRmved == T) {
  
  varNames <- openGeno[,2]
  openGeno <- openGeno[,-c(1,2)]

  nbOfSitesPerCell <- unlist(lapply(1:ncol(openGeno), function(x) { sum(!is.na(openGeno[,x])) }))
  #hist(nbOfSitesPerCell)
  quant95 <- quantile(nbOfSitesPerCell, probs=c(0.025, 0.85))
  #hist(nbOfSitesPerCell[which(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2])])
  rmCells <- which(!(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2]))
  openGeno <- openGeno[,-rmCells]

} else {
  
  varNames <- openGeno[,2]
  openGeno <- openGeno[,-c(1,2)]

}

subSet_SC <- openGeno[match(hapFile, varNames),]
colnames(subSet_SC) <- colnames(openGeno)
dim(subSet_SC)

subVarNames <- varNames[match(hapFile, varNames)]


nrow(subSet_SC)
length(hapFile)

#windowsAcrossCells <- list()
#propAcrossCells <- list()

pdf(file=outfileName, height = 11, width = 6)
par(mfrow=c(5,1))
cellNb <- 1

while (cellNb <= ncol(subSet_SC)) {
  
  index_to_compare <- which(!is.na(subSet_SC[,cellNb]))
  window <- 1
  
  
  hapOne <- c()
  hapTwo <- c()
  listWindows <- list()
  
  while (window+9 <= length(index_to_compare)) {
    
    vectPerCell <- subSet_SC[index_to_compare[window:(window+9)],cellNb]
    namesWindow <- subVarNames[index_to_compare[window:(window+9)]]
    
    names(vectPerCell) <- namesWindow
    
    subHapList <- openHap[,match(namesWindow, colnames(openHap))]
    
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
