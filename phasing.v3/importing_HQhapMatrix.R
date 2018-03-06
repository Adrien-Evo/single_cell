pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/HQhap_matrix"
fileNames <- list.files(pathToScratch, pattern="*.txt", full.names=TRUE)

openMatrix <- read.table(fileNames[1], header = T)
namesOpenMatrix <- scan(fileNames[1], what = numeric(), nlines = 1)


creatingHapList <- function(HQm) { 

  hp_l <- list()

  for (line in 1:nrow(HQm)) {
    
    newLine <- HQm[line,][!is.na(HQm[line,])]
    namesNewLine <-  colnames(HQm)[!is.na(HQm[line,])]
    
    names(newLine) <- namesNewLine
    
    hp_l[[line]] <- newLine
    
  }

  return(hp_l)
}

