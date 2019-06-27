#!/usr/bin/env Rscript

#pathToScratch <- getwd()
#pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/AD393/scriptsPhasing"
#including functions
pathToScratch <- "/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"
.libPaths( c("/.mounts/labs/awadallalab/private/flamaze/R_packages", .libPaths() ) )
library(HMM)

args <- commandArgs(trailingOnly = TRUE)

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

#################
#Function
#################
colPolygon <- function(coordTmp) {
  
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
##------
#----


folderName <- paste(pathToScratch, "/",indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, "_", dateToday, sep = "")
fileName <- paste(indId, ".", chr, ".finalHapNeighborHaps_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, ".txt", sep = "")
outfileName_M <- paste(folderName, "/InferredHapsHMap_", indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, ".txt", sep = "")
plotOutfileName <- paste(folderName, "/InferredHapsHMap_", indId, "_", chr, "_QCvar", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks, ".pdf", sep = "")
#opening haplotype file
hapFile <- scan(paste(folderName, "/", fileName, sep = ""), nlines=1, what=numeric())
openHap <- read.table(paste(folderName, "/", fileName, sep = ""), header=T, sep="\t")
colnames(openHap) <- hapFile

#GENO MATRIX FILE NAMES
prefox <- paste0("/.mounts/labs/awadallalab/scratch/SingleCell/vcfs_sc/allSamples", "/", indId, "/", indId, ".SingleCellsHetSNPs.", chr)
#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
openGeno <- data.frame(read.table(paste0(prefox, ".GT.FORMAT"), header = T, na.strings = "."))
openGenoq <- data.frame(read.table(paste0(prefox, ".GQ.FORMAT"), header = T, na.strings = "."))

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
#dim(subSet_SC)
subVarNames <- varNames[match(hapFile, varNames)]
# nrow(subSet_SC)
# length(hapFile)

#####################
##
## calling HMM
##
#####################
library(HMM)

states <- c("P", "M", "U")
symbols <- c("P1", "P2")
state_trans_prob <- matrix(c(0.9999, 0.00005, 0.00005, 0.00005, 0.9999, 0.00005, 0.00005, 0.00005,0.9999),3)
symbols_emiss_prob <- matrix(c(0.85, 0.15, 0.15, 0.85, 0.5, 0.5),nrow = 3, byrow = T)

hmm = initHMM(states, symbols, transProbs=state_trans_prob,
              emissionProbs=symbols_emiss_prob)
print(hmm)

inferredHap <- list()
index_to_compare <- list()
finalInference <- list()

for (cellNb in 1:ncol(subSet_SC)) {
  
  index_to_compare[[cellNb]] <- which(!is.na(subSet_SC[,cellNb]))
  finalInference[[cellNb]] <- rep(".", nrow(subSet_SC))
  
  subsetting_geno_matrix <- subSet_SC[index_to_compare[[cellNb]],cellNb]
  subsetting_var_names <- subVarNames[index_to_compare[[cellNb]]]
  
  
  subHapList <- openHap[,match(subsetting_var_names, hapFile)]
  
  transf_SC_hap <- subsetting_geno_matrix
  
  transf_SC_hap[unlist(lapply(1:length(subsetting_var_names), function(x) { subsetting_geno_matrix[x] == subHapList[1,x]  }))] <- "P1"
  transf_SC_hap[unlist(lapply(1:length(subsetting_var_names), function(x) { subsetting_geno_matrix[x] == subHapList[2,x]  }))] <- "P2"
  
  
  inferredHap[[cellNb]] <- viterbi(hmm, observation = transf_SC_hap)
  finalInference[[cellNb]][index_to_compare[[cellNb]]] <- inferredHap[[cellNb]]
}

df_haps <- data.frame(matrix(unlist(finalInference), ncol = ncol(subSet_SC), byrow = F))
df_haps_final <- cbind(rep(chr, length(subVarNames)), subVarNames, df_haps)
colnames(df_haps_final) <- c("Chrm", "Position", paste("cell", 1:ncol(subSet_SC), sep = "_"))
write.table(df_haps_final, file=outfileName_M, quote = F, row.names = F, col.names = T, sep="\t")

#####################                                                                                                                                      
##
## instanciate maternal and paternal cells                                                                                                                                  
##
#####################  

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


#####################
##
## drawing blocks
##
#####################
pdf(file=plotOutfileName, height = 10, width = 8)
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
  lastCoord <- which(v == "P" | v == "M" | v == "U")[length(which(v == "P" | v == "M" | v == "U"))]
  x_coord <- sort(c(mCoord, pCoord, uCoord, lastCoord))
  nbPolygons[[cellNb]] <- length(x_coord)-1
  
  
  for (polNb in 1:nbPolygons[[cellNb]]) {
    
    polygon(x=c(x_coord[polNb],x_coord[polNb],x_coord[polNb+1],x_coord[polNb+1]), y=c((countSC-1),countSC,countSC,(countSC-1)), col=colPolygon(x_coord[polNb]), border = NA)  
  }
  countSC <- countSC+1
}

dev.off()
##------------
#---------


