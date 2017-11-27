#including functions
#pathToScratch <- "/.mounts/labs/awadallalab/scratch/ialves/scriptsPhasing"
pathToScratch <- getwd()

source(paste(pathToScratch,"/building_final_haps_functions.R", sep=""))
source(paste(pathToScratch,"/computingLinks_functions.R", sep=""))
source(paste(pathToScratch,"/subsettingLinksMatrix_function.R", sep=""))
source(paste(pathToScratch,"/list_to_matrix_to_list_function.R", sep=""))
source(paste(pathToScratch,"/creatingHapList_gettingHapNames_functions.R", sep=""))
source(paste(pathToScratch,"/excludingDuplicates_function.R", sep=""))
source(paste(pathToScratch,"/mergingHap_function.R", sep=""))
source(paste(pathToScratch,"/building_final_haps_functions.R", sep=""))
source(paste(pathToScratch,"/creatingLQ_HQhap_tags_function.R", sep=""))
# .libPaths(c("/.mounts/labs/awadallalab/private/flamaze/R_packages", .libPaths()))
# library(Rmpi)
library(parallel)
dateToday <- "aug22"


#chromosome
chr <- "chr15"
#ind ID
indId <- "AD393"
#HQ ratio difference to accept link
HQ_ratio <- 5
#LQ ratio difference to accept link
LQ_ratio <- 3
#export HQ haplotypes
exportMatPhaseONe <- T
#export HQ links table
exportHQLinksTbl <- F
#export LQ links table
exportLQLinksTbl <- T
#export HQ link table
printHQcountM <- T
#nb of LQ links supporting a HQ_HQ connection
minLQsupport <- 3
##-----------
###################
# Quality filter parameters
###################

minHQCells <- 18
minNbCells <- 16
minNbLinks <- 8
##----------


#FILE NAMES
prefox <- paste(pathToScratch, "/", indId, ".genotypeMatrix.", chr, ".", sep="")
#list of heterozygous SNPs overlapping between genotyping and SC vcf file
hetSNPsList <- paste(pathToScratch, "/", indId, ".hetCallsHQ.", chr, ".c.gt", sep="")
#computation output
output <- paste(pathToScratch, "/",indId, ".", chr, ".phasing_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".log", sep="")
cat("", file=output, sep="")
#output Matrix for further analysis
linksOutput <- paste(pathToScratch, "/",indId, ".",chr, ".phaseOneHQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".txt", sep="")
#cat("", file=linksOutput, sep = "")
#output matrix with all the haplotypes
fNameMatrixSortedHap <- paste(pathToScratch, "/",indId, ".", chr, ".HQhapMatrix", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".txt", sep="")
#output all pairwise combinations phase two
phaseTwoLinksOutput <- paste(pathToScratch, "/",indId, ".",chr, ".phaseTwoHQ", minHQCells, "_NbCells", minNbCells, "_NbLinks", minNbLinks,".txt", sep="")
#log file for phasing LQ variants
outputLQPhasing <- paste(pathToScratch, "/",indId, ".",chr, ".finalPhasing_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".log", sep="")
cat("", file=outputLQPhasing, sep = "")
#log file for bridging HQs
HQ_HQ_connectionLog <- paste(pathToScratch, "/",indId, ".",chr, ".HQ_HQconnection", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".log", sep="")
cat("", file=HQ_HQ_connectionLog, sep = "")
#log file for link computation with LQ
LQlinksLog <- paste(pathToScratch, "/",indId, ".",chr, ".LQlinks_log", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks,".log", sep="")
cat("", file=LQlinksLog, sep = "")

##-----------

#opening geno and genoq files generated from subsetting the SC vcf file by extracting hetSNPsList of sites
openGeno <- data.frame(read.table(paste(prefox, "geno", sep=""), header = F, na.strings = "."))
openGenoq <- data.frame(read.table(paste(prefox, "genoq", sep=""), header = F, na.strings = "."))
openHetSNPs <- data.frame(read.table(hetSNPsList, header=T))
# namesHapPos <- scan(file=fNameMatrixSortedHap, what = numeric(), nlines = 1)
# matrixSortedHap <- data.frame(read.table(fNameMatrixSortedHap, header=T))
# names(matrixSortedHap) <- namesHapPos

#keeping positions column
hetSNPsPos <- openHetSNPs[,2]
varNames <- openGeno[,2]
openGeno <- openGeno[,-c(1,2)]
openGenoq <- openGenoq[,-c(1,2)]
nbOfCellCovPerSite <- unlist(lapply(1:nrow(openGeno), function(x) { sum(!is.na(openGeno[x,])) }))
#hist(nbOfCellCovPerSite)
nbOfSitesPerCell <- unlist(lapply(3:ncol(openGeno), function(x) { sum(!is.na(openGeno[,x])) }))
#hist(nbOfSitesPerCell)
quant95 <- quantile(nbOfSitesPerCell, probs=c(0.025, 0.975))
#hist(nbOfSitesPerCell[which(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2])])
rmCells <- which(!(nbOfSitesPerCell > quant95[1] & nbOfSitesPerCell < quant95[2]))
openGeno <- openGeno[,-rmCells]
openGenoq <- openGenoq[,-rmCells]

# openGeno <- openGeno[,-c(1,2)]
# openGenoq <- openGenoq[,-c(1,2)]

#######################
##
## Compute all pairwise comb for PHASE ONE 
##
#######################


lengthHQ <- apply(openGenoq, 1, function(x) { length(which(x > 20))})
HQvarNames <- varNames[which(lengthHQ >= minHQCells)] #replace 25 by minHQCells
HQ_m <- openGeno[which(lengthHQ >= minHQCells),] #replace 25 by minHQCells
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#cl <- makeCluster(no_cores, type="FORK")
# Initiate cluster
#cl <- makeCluster( mpi.universe.size(), type = 'MPI')
clusterExport(cl, c("minNbCells", "get_LinkCount", "computingLinks", "chr", "HQvarNames", "HQ_m", "output"))
listLinks <- list()
#creating a file with the links and corresponding counts
#this file "linksOutput" is required for step 2 below
#computingLinks <- function(tmpSite, listOfVar, HQ_genotypeMatrix, outFileName) 
listLinks <- parLapply(cl=cl, HQvarNames, function(x) { computingLinks(tmpSite=x, listOfVar=HQvarNames, HQ_genotypeMatrix=HQ_m, outFileName=output) })
stopCluster(cl)
fullHQlinksComb <- data.frame()
for (HQvar in 1:length(listLinks)) {
  
  if (length(listLinks[[HQvar]]) > 1) {
    
    if (HQvar == 1) {
      fullHQlinksComb <- listLinks[[HQvar]]
    } else {
      fullHQlinksComb <- rbind(fullHQlinksComb,listLinks[[HQvar]])
    }
    
  } 
}

links_and_order_tmp <- subsettingLinksMatrix(fullHQlinksComb, HQ_ratio)
subMatrix <- links_and_order_tmp[[1]]
subOrder <- links_and_order_tmp[[2]]

if (printHQcountM) {
  
  write.table(subMatrix, file=linksOutput, quote = F, row.names = F, col.names = T, sep="\t", append = T)
  
}


rm(links_and_order_tmp)
rm(fullHQlinksComb)
#creating a matrix with all haplotypes 
m_hap <- creatingHapMatrixFromLinkCountMatrix(subMatrix,subOrder)

#transforming m_hap into a list
l_hap <- creatingHapListFromLinkCountMatrix(m_hap)

#### COMMENT OUT if you are generating the hapList
# scanNames <- scan(fNameMatrixSortedHap, nlines = 1,what = character(), sep = "\t")
# tmpM <- read.table(fNameMatrixSortedHap, header=T)
# colnames(tmpM) <- scanNames

l_hap <- creatingHapListFromLinkCountMatrix(tmpM)
###----------
#---

#the next steps remove direct duplicates followed by merging them (and merging partial duplicates - sharing the same var Names but with the inverted
#haplotypes)
uniqueHap <- lapply(1:length(l_hap), function(x) { creatingHapList(x, l_hap) })
uniqueHap <- unique(uniqueHap)
#generating a list of unique merged haplotypes: FINAL config
tmpListToRm <- excludingDuplicates(uniqueHap)
tmpListToRm <- lapply(1:length(tmpListToRm), function(x) { creatingHapList(x, tmpListToRm) })
#rm(uniqueHap)


hapList <- tmpListToRm
namesHapList <- namesHapListFunction(hapList)
rm(tmpListToRm)
rm(subMatrix)
rm(subOrder)


#convert all elements of the hapList into the same type of objects
hapList <- lapply(1:length(hapList), function(x) { creatingHapList(x, hapList) })

if (exportMatPhaseONe == T) {
  
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


#resetting QC filters
minHQCells <- 10
minNbCells <- 4
minNbLinks <- 2 


varInHQhaps <- namesHapList
rm(namesHapList)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#cl <- makeCluster(no_cores, type="FORK")
# Initiate cluster
#cl <- makeCluster( mpi.universe.size(), type = 'MPI')
clusterExport(cl, c("minNbCells", "get_LinkCount", "computingLinks", "chr", "openGenoq", "openGeno", "varNames", "outputLQPhasing", "minHQCells", "varInHQhaps", "subsettingByQuality_computingLQLinks"))
LQ_HQ_links_list <- list()

#openGenoq, openGeno, varNames, outputLQPhasing, minHQCells

LQ_HQ_links_list <- parLapply(cl=cl, varNames, function(x) { subsettingByQuality_computingLQLinks(tmpLQSite=x, varNamesInHQhaps=varInHQhaps, outLog=LQlinksLog) } )
#change from varNames to varNames-HQvar
stopCluster(cl)


fullHQlinksComb <- data.frame()
for (HQvar in 1:length(LQ_HQ_links_list)) {
  
  if (length(LQ_HQ_links_list[[HQvar]]) > 1) {
    
    if (HQvar == 1) {
      fullHQlinksComb <- LQ_HQ_links_list[[HQvar]]
    } else {
      fullHQlinksComb <- rbind(fullHQlinksComb,LQ_HQ_links_list[[HQvar]])
    }
    
  } 
}

#### COMMENT OUT if you are generating the phaseTwoLinksOutput matrix
#fullHQlinksComb <- read.table(phaseTwoLinksOutput, header=T)

###--------
#-----

links_and_order_tmp <- subsettingLinksMatrix(fullHQlinksComb, LQ_ratio)
subMatrix <- links_and_order_tmp[[1]]
subOrder <- links_and_order_tmp[[2]]
rm(links_and_order_tmp)
rm(fullHQlinksComb)

varLQstep <- names(table(subMatrix[,2]))


if (exportLQLinksTbl) {
  
  write.table(subMatrix, file=phaseTwoLinksOutput, quote = F, row.names = F, col.names = T, sep="\t")
  
}

## Calculate the number of cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#cl <- makeCluster(no_cores, type="FORK")
# Initiate cluster
#cl <- makeCluster( mpi.universe.size(), type = 'MPI')
clusterExport(cl, c("hapList", "subMatrix", "subOrder","outputLQPhasing","creatingLQhapsWithTags", "taggingLQvar_HQhap_combination", "varLQstep"))
finalHapList <- list()

finalHapList <- parLapply(cl=cl, 1:length(varLQstep), function(varPos) { creatingLQhapsWithTags(varLQstep[varPos]) })
 
stopCluster(cl)   


allHapsTagsTogether <- unlist(finalHapList)
allHapsTagsTogether <- allHapsTagsTogether[which(allHapsTagsTogether != 0)]
hap_df <- data.frame(matrix(unlist(strsplit(allHapsTagsTogether, split = "_")), ncol = 5, byrow = T))
#HQlinks_support <- lapply(1:(length(hapList)-1), function(l) { hap_df[which(hap_df[,2] == l & hap_df[,3] == l+1),] })
connectionTag <- buildingHQ_connections(length(hapList),outF=HQ_HQ_connectionLog, t=minLQsupport)
listHQhap_LQvar <- linkingLQvar_toHQhaplotypes()

perHQhapCombinations <- lapply(1:length(listHQhap_LQvar), function(z) assemblingCompleteHaplotypes(z))


hapListID_toMerge <- which(!is.element(1:length(hapList), as.numeric(unique(as.vector(matrix(unlist(strsplit(connectionTag, split = "_")), ncol=3, byrow = T)[,1:2])))))
if (length(hapListID_toMerge) > 0) {
  perHQhapCombinations[c(hapListID_toMerge)] <- 0 
}  

fullHapZero <- c()
fullHapOne <- c()
firstTag <- as.numeric(unlist(strsplit(connectionTag[1], split="_"))[1])

fullHapZero <- perHQhapCombinations[[firstTag]][[1]]
fullHapOne <- perHQhapCombinations[[firstTag]][[2]]
HQcomb <- firstTag+1


while(HQcomb <= length(connectionTag)+1) {

#HQcomb <- 2
tmpTag <- unlist(strsplit(connectionTag[HQcomb-1], split="_"))[3]
hapIndex <- unlist(strsplit(connectionTag[HQcomb-1], split="_"))[1:2]
print(paste("Including subhaplotype nb:", hapIndex[2], sep = " "))
      
    if (tmpTag == "AD") {
      
      fullHapZero <- c(fullHapZero, perHQhapCombinations[[as.numeric(hapIndex[2])]][[1]])
      fullHapOne <- c(fullHapOne, perHQhapCombinations[[as.numeric(hapIndex[2])]][[2]])
      
    } else if (tmpTag == "BC") {
      
      fullHapZero <- c(fullHapZero, perHQhapCombinations[[as.numeric(hapIndex[2])]][[2]])
      fullHapOne <- c(fullHapOne, perHQhapCombinations[[as.numeric(hapIndex[2])]][[1]])
      
    }

HQcomb <- HQcomb+1
}
fullHaplotypes_df <- rbind(fullHapZero,fullHapOne)
orderedFullHaplotypes_df <- fullHaplotypes_df[,order(as.numeric(colnames(fullHaplotypes_df)), decreasing = F)]

write.table(orderedFullHaplotypes_df, file=paste("fullHaplotypes_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".txt", sep=""), quote = F, row.names = F, col.names = T, sep = "\t")





#######################
##
##
#######################