#!/usr/bin/env Rscript


##########################
#Functions include breakpoint_computation, size_unknown, size_M, Size_P -> stored in R script breakpoint_functions.R
#Global variable -> vcf file that is stored as "df"
#Functions take df as an argument, which specifies the indiviudal and chromosome it orginated from
#Calculates the sizes of the breakpoint, unknown, M and P blocks in every cell of the file (cell = column)
#Outputs the results into four corresponding histograms
#output filename (.pdf) = input filename (.txt) with different extension

#Created by Charlotte Michelle Nguyen
#November, 2017

##########################


pathToFile <- "/home/isabelalves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/Oct8/resolution/"
#pathToFile <- getwd()
#import library to execute functions in parallel 
library(parallel)
source(paste(pathToFile,"/breakpoint_functions.R", sep = ""))

#ADD commands defining the data frame file name
# open the data frame and save is as df
dataframes <- list.files(pattern = "*chr21_QCvar25_NbCells10_NbLinks5.txt$")
bp_size_list <- list()
nb_of_breakpoints <- list()
nb_of_cells_no_bp <- list()
nb_of_cells <- list()
medSizeBreakpoint <- list()
bp_size_list_clean <- list()
medSizeBreakpoint_clean <- list()
q.lim <- .95
#ADD output file specifications



for (i in 1:length(dataframes)) {
  
  i <- 1
  df <- read.table(file = dataframes[i], header = TRUE, stringsAsFactors = FALSE)
  outputFile <- sub("txt", "pdf", dataframes[i])
  indID <- unlist(strsplit(unlist(strsplit(dataframes[i], split="\\."))[1], split = "_"))[2]
  chrID <- unlist(strsplit(unlist(strsplit(dataframes[i], split="\\."))[1], split = "_"))[3]
  print(paste0(indID, "_", chrID))
  #Takes vcf file specifying individual and chromosome # to detect breakpoints within single cells 
  #df <- AD393.chr15.inferredHaps_HQ25_NbCells20_NbLinks10_RmCells
  
  #List of breakpoints, unknown, M, and P blocks
  list_breakpoint_size <- list()
  list_unknown_size <- list()
  list_M_size <- list()
  list_P_size <- list()
  
  #parLapply -> speeds up functions by splitting the computation in paralell 
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  clusterExport(cl, c("df", "breakpoint_computation", "size_unknown", "size_M", "size_P"))
  list_breakpoint_size <- parLapply(cl, 3:ncol(df), function(x) { breakpoint_computation(df[,x]) })
  list_unknown_size <- parLapply(cl, 3:ncol(df), function(x) { size_unknown(df[,x]) })
  list_M_size <- parLapply(cl, 3:ncol(df), function(x) { size_M(df[,x]) })
  list_P_size <- parLapply(cl, 3:ncol(df), function(x) { size_P(df[,x]) })
  stopCluster(cl)
  
  nb_of_breakpoints[[i]] <- lengths(list_breakpoint_size)
  bp_size_list[[i]] <- unlist(list_breakpoint_size) #the  NULL entries disapear
  
  nb_of_cells_no_bp[[i]] <- sum(unlist(lapply(list_breakpoint_size, is.null))) #changed by Dec 3
  nb_of_cells[[i]] <- length(list_breakpoint_size) #this takes into account the NULL entries
  medSizeBreakpoint[[i]]  <- median(bp_size_list[[i]][bp_size_list[[i]] > 0]/1000)

 
  bp_size_list_clean[[i]] <- unlist(list_breakpoint_size[-c(which(nb_of_breakpoints[[i]] > quantile(nb_of_breakpoints[[i]], probs = q.lim)))])
  medSizeBreakpoint_clean[[i]] <-  median(bp_size_list[[i]][-c(which(nb_of_breakpoints[[i]] > quantile(nb_of_breakpoints[[i]], probs = q.lim)))]/1000)

  list_breakpoint_size_clean <- list_breakpoint_size[-c(which(nb_of_breakpoints[[i]] > quantile(nb_of_breakpoints[[i]], probs = q.lim)))]

  #Applies functions to every cell in vcf file, stores output in different list
  #list_breakpoint_size <- lapply(3:ncol(df), function(x) { breakpoint_computation(df[,x]) })
  #list_unknown_size <- lapply(3:ncol(df), function(x) { size_unknown(df[,x]) })
  #list_M_size <- lapply(3:ncol(df), function(x) { size_M(df[,x]) })
  #list_P_size <- lapply(3:ncol(df), function(x) { size_P(df[,x]) })
   #Make histogram of size of breakpoints
  histDist <- hist(bp_size_list_clean[[i]][bp_size_list_clean[[i]] > 0]/1000,
       main="Resolution of Breakpoints",
       xlab="Size of interval (kb)",
       ylab = "Number of breakpoints", breaks = c(seq(0, max(bp_size_list_clean[[i]])/1000, by = 100), max(bp_size_list_clean[[i]])+10/1000), plot = F
  )
  cdist <- ecdf(bp_size_list_clean[[i]][bp_size_list_clean[[i]] > 0]/1000)
  barplot(histDist$counts,  
          xlab="Size of interval (kb)",
        ylab = "Number of breakpoints", xaxt='n', space = 0, main = paste0(indID, "_", chrID))
  lines(x = c(0:length(histDist$counts)), y=c(0,cdist(histDist$mids)*max(histDist$counts)), col ='red')
  axis(4, at=seq(from = 0, to = max(histDist$counts), length.out = 11), labels=seq(0, 1, 0.1), col = 'red', col.axis = 'red')
  axis(1, at=seq(0,length(histDist$counts), length.out = length(c(seq(0, max(bp_size_list_clean[[i]])/1000, by = 500)))+1), 
                 labels=c(seq(0, max(bp_size_list_clean[[i]])/1000, by = 500), (max(bp_size_list_clean[[i]])+10)/1000), col = 'black', col.axis = 'black')
  mtext(side = 4, line = 3, 'Cumulative Density', col = 'red')
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts), legend=paste0("Mean nb of breakpoints: ", round(sum(nb_of_breakpoints[[i]][-c(which(nb_of_breakpoints[[i]] > quantile(nb_of_breakpoints[[i]], probs = q.lim)))])/nb_of_cells[[i]], digits = 2)), bty = "n")  
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts)-(max(histDist$counts)/6), legend=paste0("Nb of cells: ", nb_of_cells[[i]]), bty = "n")  
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts)-(max(histDist$counts)/4), legend=paste0("Median size bp: ", round(medSizeBreakpoint_clean[[i]], digits = 2)), bty = "n")  
  dev.off()

}


