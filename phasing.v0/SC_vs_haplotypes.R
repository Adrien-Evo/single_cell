pathToScratch <- getwd()
pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/scriptsPhasing"
#opening haplotype file
hapFile <- scan("fullHaplotypes_Aug22_HQ18_NbCells16_NbLinks8.txt", nlines=1, what=numeric())
openHap <- read.table(file="fullHaplotypes_Aug22_HQ18_NbCells16_NbLinks8.txt", header=T, sep="\t")
colnames(openHap) <- hapFile

indId <- "AD393"
chr <- "chr15"
popPhased <- FALSE
if (popPhased) {
  
  m_phased_name <- "C:/Users/ialves/Dropbox/singleCellProject/phasingConfirmation_GT/Omni_SNPchip__chr15.ma"
  openM_phased <- read.table(m_phased_name, header=T, sep="\t")
  
  hapOne <- unlist(lapply(1:nrow(openM_phased), function(x) {unlist(strsplit(as.character(openM_phased[x,10]), split="|"))[1]}))
  haptwo <- unlist(lapply(1:nrow(openM_phased), function(x) {unlist(strsplit(as.character(openM_phased[x,10]), split="|"))[3]}))
  
  openHap <- rbind(as.numeric(hapOne), as.numeric(haptwo))
  colnames(openHap) <- openM_phased[,2]
  hapFile <- colnames(openHap)
  
}
#FILE NAMES
prefox <- paste(pathToScratch, "/", indId, ".genotypeMatrix.", chr, ".", sep="")
#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
openGeno <- data.frame(read.table(paste(prefox, "geno", sep=""), header = F, na.strings = "."))
openGenoq <- data.frame(read.table(paste(prefox, "genoq", sep=""), header = F, na.strings = "."))

varNames <- openGeno[,2]
openGeno <- openGeno[,-c(1,2)]
openGenoq <- openGenoq[,-c(1,2)]

nbOfSitesPerCell <- unlist(lapply(1:ncol(openGeno), function(x) { sum(!is.na(openGeno[,x])) }))
hist(nbOfSitesPerCell)
quant95 <- quantile(nbOfSitesPerCell, probs=c(0.025, 0.975))
hist(nbOfSitesPerCell[which(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2])])
rmCells <- which(!(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2]))
openGeno <- openGeno[,-rmCells]



subSet_SC <- openGeno[match(hapFile, varNames),]
colnames(subSet_SC) <- colnames(openGeno)
dim(subSet_SC)

subVarNames <- varNames[match(hapFile, varNames)]


nrow(subSet_SC)
length(hapFile)

windowsAcrossCells <- list()
propAcrossCells <- list()

pdf(file="SCcomp_phasedHap_Aug22.pdf", height = 11, width = 6)
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
  windowsAcrossCells[[cellNb]] <- listWindows
  propAcrossCells[[cellNb]] <- rbind(hapOne,hapTwo) 
  cellNb <- cellNb+1
}

dev.off()

#getting indexes of sites with information from the SC

#comparing with pop phased data
#openHap <- m_popPhased
# 
# unlist(lapply(1:length(subsetting_var_names), function(x) { subsetting_geno_matrix[x] == subHapList[1,x]  }))
# unlist(lapply(1:length(subsetting_var_names), function(x) { subsetting_geno_matrix[x] == subHapList[2,x]  }))

install.packages("HMM")
library(HMM)

states <- c("P", "M", "U")
symbols <- c("P1", "P2")
state_trans_prob <- matrix(c(0.9999, 0.00005, 0.00005, 0.00005, 0.9999, 0.00005, 0.00005, 0.00005,0.9999),3)
symbols_emiss_prob <- matrix(c(0.90, 0.10, 0.10, 0.90, 0.5, 0.5),nrow = 3, byrow = T)

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
df_haps_final <- cbind(rep("chr15", length(subVarNames)), subVarNames, df_haps)
colnames(df_haps_final) <- c("Chrm", "Position", paste("cell", 1:ncol(subSet_SC), sep = "_"))
write.table(df_haps_final, file="inferredHaplotypes_H20_Aug18.txt", quote = F, row.names = F, col.names = T, sep="\t")

