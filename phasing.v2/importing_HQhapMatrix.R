pathToScratch <- "C:/Users/ialves/Dropbox/singleCellProject/phasing_donors/HQhap_matrix"
fileNames <- list.files(pathToScratch, pattern="*.txt", full.names=TRUE)

openMatrix <- read.table(fileNames[1], header = T)
namesOpenMatrix <- scan(fileNames[1], what = numeric(), nlines = 1)
hapList <- list()

for (line in 1:nrow(openMatrix)) {
  
  newLine <- openMatrix[line,][!is.na(openMatrix[line,])]
  namesNewLine <-  namesOpenMatrix[!is.na(openMatrix[line,])]
  
  names(newLine) <- namesNewLine
  
  hapList[[line]] <- newLine
  
}


