
hapFile <- scan("fullHaplotypes_Aug22_HQ18_NbCells16_NbLinks8.txt", nlines=1, what=numeric())
openHap <- read.table(file="fullHaplotypes_Aug22_HQ18_NbCells16_NbLinks8.txt", header=T, sep="\t")
colnames(openHap) <- hapFile



###################
# Quality filter parameters
###################
minHQCells <- 8
minNbCells <- 4
minNbLinks <- 2


randomVar <- sort(sample(hapFile, 100), decreasing = F)

compareWithVar <- sort(sample(hapFile[which(hapFile != randomVar[2])], 100), decreasing = F)


toSubSet <- c(hapFile[1], hapFile[2:10])

HQvarNames <- varNames[match(toSubSet, varNames)] #replace 25 by minHQCells
HQ_m <- openGeno[match(toSubSet, varNames),] 
outF_error <- paste(pathToScratch, "errorRates_HQ", minHQCells, "_NbCells", minNbCells, "_NbLinks",minNbLinks, ".txt", sep = "")


computingLinks(hapFile[1],HQvarNames,HQ_m,outF_error)



