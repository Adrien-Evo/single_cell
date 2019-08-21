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


pathToFile <- "/.mounts/labs/awadallalab/public/ialves/HMM_Jan29/"
pathToFile <- getwd()
#import library to execute functions in parallel 
library(parallel)
source(paste(pathToFile,"/breakpoint_functions.R", sep = ""))

#ADD commands defining the data frame file name
# open the data frame and save is as df
dataframes <- list.files(pattern = "*.txt$")
bp_size_list <- list()

for (i in 1:length(dataframes)) {
  
  i <- 3
  df <- read.table(file = dataframes[i], header = TRUE, stringsAsFactors = FALSE)
  outputFile <- sub("txt", "pdf", dataframes[i])


    
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
  
  bp_size_list[[i]] <- unlist(list_breakpoint_size)
  #Applies functions to every cell in vcf file, stores output in different list
  #list_breakpoint_size <- lapply(3:ncol(df), function(x) { breakpoint_computation(df[,x]) })
  #list_unknown_size <- lapply(3:ncol(df), function(x) { size_unknown(df[,x]) })
  #list_M_size <- lapply(3:ncol(df), function(x) { size_M(df[,x]) })
  #list_P_size <- lapply(3:ncol(df), function(x) { size_P(df[,x]) })
  
  #ADD output file specifications
  pdf(file = outputFile, height = , width = )
  par(mfrow=c(2,2))
  
  #Make histogram of size of breakpoints
  hist(unlist(list_breakpoint_size)[unlist(list_breakpoint_size) > 0]/1000,
       main="Resolution of Breakpoints",
       xlab="Size of interval (kb)",
       ylab = "Number of breakpoints"
  )
  
  #Make histogram of size of unknown blocks
  hist(unlist(list_unknown_size)[unlist(list_unknown_size) > 0]/1000,
       main="Resolution of Unknown Blocks",
       xlab="Size of interval (kb)",
       ylab = "Number of unknown blocks"
  )
  
  #Make histogram of size of M blocks
  hist(unlist(list_M_size)[unlist(list_M_size) > 0]/1000,
       main="Resolution of M Blocks",
       xlab="Size of interval (kb)",
       ylab = "Number of M blocks"
  )
  
  #Make histogram of size of P blocks
  hist(unlist(list_P_size)[unlist(list_P_size) > 0]/1000,
       main="Resolution of P Blocks",
       xlab="Size of interval (kb)",
       ylab = "Number of P blocks"
  )
  
  
  dev.off()
}
