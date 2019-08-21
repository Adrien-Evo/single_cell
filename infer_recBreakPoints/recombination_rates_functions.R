#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

##########################
#Functions include breakpoint_chrom_position, unknown_chrom_position, breakpoint_counts, cell_counts -> stored in R script recombination_rates_functions.R
#Global variable -> vcf file that is stored as "df"
#Functions take df as an argument, which specifies the indiviudal and chromosome it orginated from
#Calculates the coordinates of breakpoints and unknown blocks, and then determines the breakpoints and cell counts per 1Mb bin


#Created by Charlotte Michelle Nguyen
#November, 2017

##########################

##########################
###
###
###    Functions
###
###
#########################
#Computes coordinates of breakpoints in every cell of file
breakpoint_chrom_position <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,109]
  
  #Stores first and last coordinate of each M/P block
  bp <- c()
  #Neutral parental state of the chromosome (vs. M/P)
  state <- "I"
  
  #Scans cell from the first to last position of chromosome to detect parental origin (M or P)
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of M block
    if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores first coordinate of M block
      bp <- c(bp, i)
      #Scans for previous P block end coordinate up until new M block
      bp <- c(bp, suppressWarnings(max(which(v[1:i] == "P"))))
      
      #Scans for first coordinate of new P block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores the first coordinate of P block 
      bp <- c(bp, i)
      #Scans for previous M block end coordinate up until new P block
      bp <- c(bp, suppressWarnings(max(which(v[1:i] == "M"))))
    }
  } #end of going along every element of the vector with the parental states
  
  #Removes all NA/null values of bp
  bp <- bp[is.finite(bp)]
  
  #Breakpoint computation requires min. 2 values to compute coordinates of breakpoint between M and P
  if(length(bp) > 1) {
    
    #Sorts coordinates in ascending order
    bp <- sort(c(bp)) 
    
    #Breakpoint requires more than 1 parental state in cell 
  } else {
    bp <- 0
  }
  
  #List of chromosome positions for breakpoint coordinates
  bpCoord <- c()
  
  #Removes first coordinate of bp if length >2 because first coordinate is the beginning of an M or P block
  #Stores cooresponding bp coordinate as chromosome position
  if (length(bp) > 2 ) {
    bpCoord <- c(bpCoord, df[,2][bp][2:length(bp)])
    #if bp coordinate = 2, store cooresponding chromosome positions of both coordinates as it represents a true breakpoint
  } else if (length(bp) == 2) {
    bpCoord <- c(bpCoord, df[,2][bp])
  } else {
    bpCoord <- 0
  }
  
  #Returns the size of different breakpoints 
  #print(bp)
  return(bpCoord) 
}
##----
#--


#Function to compute size of unknown blocks
unknown_chrom_position <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,109]
  
  #uCoord stores first and last coordinates of unknown blocks
  uCoord <- c()
  #Neutral parental state of the chromosome
  state <- "I"
  
  #Scans cell from  first to last position of chromosome to detect unknown positions
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of new U block
    if (v[i] != "." & v[i] == "U" & v[i] != state) {
      state <- "U"
      #Stores first cooridnate of U block 
      uCoord <- c(uCoord, i) 
      
      #Scans for new M block and finds previous end point of U block
    } else if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores U block end coordinate
      uCoord <- c(uCoord, suppressWarnings(max(which(v[1:i] == "U"))))
      
      #Scans for new P block and finds previous end point of U block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores U block end cooridnate in uCoord
      uCoord <- c(uCoord, suppressWarnings(max(which(v[1:i] == "U"))))
      
      #If U is the last state of chromosome, finds end point of U block
    } else if (i == length(v) & state == "U") {
      uCoord <- c(uCoord, suppressWarnings(max(which(v == "U"))))
    }
    
  } #end of going along every element of the vector with the parental states
  
  
  #Removes all NA/null values of uCoord and duplicate values
  uCoord <- unique(uCoord[is.finite(uCoord)])
  uChrom <- c()
  
  #size_unknown computation requires min. 2 values to compute distance between first and last U coordinate
  if(length(uCoord) > 1) {
    
    #Sorts coordinates of uCoord in ascending order
    #uCoord <- sort(c(uCoord))
    #print(ordereduCoord)
    
    #Stores size of unknown blocks, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the first and last coordinates of U blocks  
    #uCoordBlocks <- c()
    #Compute difference between last and first point of U blocks
    #for (i in seq(1, (length(uCoord)-1), 2)) {
    #uCoordBlocks <- c(uCoordBlocks, (df[,2][uCoord[i+1]] - df[,2][uCoord[i]]))
    uChrom <- c(uChrom, df[,2][uCoord]) 
    
    #No size_unknown computation if there are no U blocks
  } else {
    
    uChrom <- 0
  }
  
  #Returns the size of different unknown blocks 
  #print(uChrom) 
  return(uChrom)
}
##---------
#-----


#Test individual cell
#cell <- list_breakpoint_coordinates[[5]]
#for (cell in list_breakpoint_coordinates) {

breakpoint_counts <- function(cell) { #cell - vector with the breakpoint coordinates within a cell; 
  #it takes the list_breakpoint_cooridnate as global variable
  #cell <- list_breakpoint_coordinates[[54]]
  #FUnction only applies to cells with breakpoints
  if (length(cell) > 1) {
    #Check each pair of breakpoint coordinates (first and last)
    for (i in seq(1, length(cell)-1, 2)) {
      #Determine bin of the first coordinate
      (max(which(bin_start < cell[i])))
      #Determine bin of the last coordinate
      (min(which(bin_end > cell [i+1])))
      #If first and last cooridnate of breakpoint are in the same bin, increment count by 1
      if ((max(which(bin_start < cell[i]))) ==  (min(which(bin_end > cell [i+1])))) {
        bin_counts[(max(which(bin_start < cell[i])))] <- bin_counts[(max(which(bin_start < cell[i])))] + 1
        #If the difference between bin_end and bin_start > 1
      } else if ( (min(which(bin_end > cell [i+1]))) - (max(which(bin_start < cell[i]))) > 1 ) {
        #Store the value of the difference between bin_end and bin_start
        x <- (min(which(bin_end > cell [i+1]))) - (max(which(bin_start < cell[i])))
        #If the difference is even, increment count of median bin
        if (x %% 2 == 0) {
          bin_counts[(max(which(bin_start < cell[i]))) + x/2] <- bin_counts[(max(which(bin_start < cell[i]))) + x/2] + 1
          #If the difference is odd, increment count of smaller ofmedian bin
        } else {
          bin_counts[(max(which(bin_start < cell[i]))) + (x-1)/2] <- bin_counts[(max(which(bin_start < cell[i]))) + x/2] + 1
        }
        #If the difference between bin_end and bin_start = 1
      } else {
        #If 50% of breakpoint is covered by first bin, increment count of first bin
        if ((cell[i+1]-cell[i])/2 + bin_start[(max(which(bin_start < cell[i])))] < bin_end[(max(which(bin_start < cell[i])))]) {
          bin_counts[(max(which(bin_start < cell[i])))] <-  bin_counts[(max(which(bin_start < cell[i])))] + 1
          #If 50% of breakpoint not covered by first bin, increment count of second bin
        } else {
          bin_counts[(min(which(bin_end > cell [i+1])))] <- bin_counts[(min(which(bin_end > cell [i+1])))] + 1
        }
      }
    }
  } else {bin_counts <- 0}
  
  
  return(bin_counts)
}

##---------
#-----

#list_tmp <- list_breakpoint_counts

#list_tmp <- list_breakpoint_counts[which(lengths(list_breakpoint_counts) > 1)]
#m_tmp <- matrix(unlist(list_tmp), ncol = length(bin_counts), byrow = T)
#apply(m_tmp, 2, sum)


#############################
###
### End of functions
###
#############################

##------------------------------

#############################
###
###         Main
###
#############################
source("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/scripts_resolution/breakpoint_functions.R")

chrNb <- as.numeric(args[1])
chrNb <- 7

print(paste0("Computing recombination break points for chrm: ", chrNb))
chrTag <- paste0("chr", chrNb)
hapMapRecombFile <- paste0("genetic_map_GRCh37_", chrTag, ".txt")
deCodeRecombFile <- paste0("male.map.", chrTag, ".txt")

if (chrNb < 8) {
  
   pathToFile_tmp  <- paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/")
   dir.create(pathToFile_tmp)
   lapply(list.files(path=paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "p/"), 
                     pattern = "InferredHapsHMap*"), function(x) { file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "p/",x),
                                                                             pathToFile_tmp, recursive = F, copy.mode = T) }) 
   lapply(list.files(path=paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "q/"), 
                     pattern = "InferredHapsHMap*"), function(x) { file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "q/",x),
                                                                             pathToFile_tmp, recursive = F, copy.mode = T) }) 
   inferredHap_tmp <- list.files(path=paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag), 
              pattern = paste0("*", chrTag, "p_QCvar25_NbCells10_NbLinks5.txt"))
   setwd(pathToFile_tmp)
   for (i in 1:length(inferredHap_tmp)) {
     indID <- unlist(strsplit(inferredHap_tmp[i], split = "_"))[2]  
     headMat_tmp1 <- scan(inferredHap_tmp[i], what = character(), nlines = 1)
     mat_tmp2_name <- paste0("InferredHapsHMap_", indID, "_", chrTag, "q_QCvar25_NbCells10_NbLinks5.txt")
     headMat_tmp2 <- scan(mat_tmp2_name, what = character(), nlines = 1)
     length(headMat_tmp1) == length(headMat_tmp2)
     if (sum(is.na(match(headMat_tmp1, headMat_tmp2))) > 0) {
       print("Chromosome arms of chr: ", chrNb, " do not have info for the same cells")
       break;
     } else {
       m1 <- read.table(file = inferredHap_tmp[i], header = TRUE, stringsAsFactors = FALSE)
       m2 <- read.table(file = mat_tmp2_name, header = TRUE, stringsAsFactors = FALSE)
       full_m_name <- paste0("InferredHapsHMap_", indID, "_", chrTag, "_QCvar25_NbCells10_NbLinks5.txt")
       write.table(rbind(m1,m2), file=full_m_name, col.names = T, row.names = F, quote = F)
       file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", hapMapRecombFile),
                 pathToFile_tmp, recursive = F, copy.mode = T)
       file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", deCodeRecombFile),
                 pathToFile_tmp, recursive = F, copy.mode = T)
       print(paste0("Full matrix of inferred chromosomes printed for chr: ", chrNb, " and ind: ", indID))
       
     }
}
}
rm(pathToFile_tmp, headMat_tmp1, headMat_tmp2, m1, m2, mat_tmp2_name)

pathToFile <- paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/resolution")
dir.create(pathToFile)
lapply(list.files(path=paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag), 
                  pattern = "InferredHapsHMap*"), function(x) { file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/",x),
                                                                          pathToFile, recursive = F, copy.mode = T) }) 
file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/",hapMapRecombFile),
          pathToFile, recursive = F, copy.mode = T)
file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/",deCodeRecombFile),
          pathToFile, recursive = F, copy.mode = T)
dir.create(paste0(pathToFile, "/", "resolution.95"))

setwd(pathToFile)
#import library to execute functions in parallel 
library(parallel)


####---------------
###-------------
#Section applies for finding breakpoints/bin 

  chrSizes <- c(249250621, 243199373,198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
                141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 
                81195210, 78077248, 59128983, 63025520, 48129895, 51304566)
  
  #Chromosome 15
  chrmLength <- chrSizes[chrNb]
  print(chrmLength)
  windowSize <- 1000000 #1Mb
  bin_start<- c(seq(1, chrmLength, by=windowSize))
  bin_end <- c(seq(windowSize, chrmLength, by=windowSize), chrmLength)
  bin_counts <- vector(mode = "numeric" , length = length(bin_start))
  nb_of_breakpoints <- list()
  list_breakpoint_counts <- list()
  q.lim <- .95
  #Chromosome 21
  # chrmLength <- chrSizes[21]
  # windowSize <- 1000000 #1Mb
  # bin_start<- c(seq(1, chrmLength, by=windowSize))
  # bin_end <- c(seq(windowSize, chrmLength, by=windowSize), chrmLength)
  # bin_counts <- vector(mode = "numeric" , length = length(bin_start))
  
  
  
  #ADD commands defining the data frame file name
  # open the data frame and save is as df
  dataframes <- list.files(pattern = paste0("*", chrTag, "_QCvar25_NbCells10_NbLinks5.txt$"))
  #for (i in 1:length(dataframes)) {
    #df <- read.table(file = dataframes[i], header = TRUE, stringsAsFactors = FALSE)
    #outputFile <- sub("txt", "pdf", dataframes[i])
  
  for (file in 1:length(dataframes)) {
  
    #Takes vcf file specifying individual and chromosome # to detect breakpoints within single cells 
    #file <- 1  
    df <- read.table(file = dataframes[file], header = TRUE, stringsAsFactors = FALSE)
    
    #List of breakpoints, unknown, M, and P blocks
    list_breakpoint_coordinates <- list()
    list_unknown_coordinates <- list()
    #List of breakpoints, unknown, M, and P blocks
    list_breakpoint_size <- list()
    list_unknown_size <- list()
    list_M_size <- list()
    list_P_size <- list()
    
    #parLapply -> speeds up functions by splitting the computation in paralell 
    no_cores <- detectCores()-1
    cl <- makeCluster(no_cores)
    clusterExport(cl, c("df", "breakpoint_chrom_position", "unknown_chrom_position"))
    list_breakpoint_coordinates <- parLapply(cl, 3:ncol(df), function(x) { breakpoint_chrom_position(df[,x]) })
    list_unknown_coordinates <- parLapply(cl, 3:ncol(df), function(x) { unknown_chrom_position(df[,x]) })
    stopCluster(cl)
    
    #}
    #parLapply -> speeds up functions by splitting the computation in paralell 
    no_cores <- detectCores()-1
    cl <- makeCluster(no_cores)
    clusterExport(cl, c("df", "breakpoint_computation", "size_unknown", "size_M", "size_P"))
    list_breakpoint_size <- parLapply(cl, 3:ncol(df), function(x) { breakpoint_computation(df[,x]) })
    list_unknown_size <- parLapply(cl, 3:ncol(df), function(x) { size_unknown(df[,x]) })
    list_M_size <- parLapply(cl, 3:ncol(df), function(x) { size_M(df[,x]) })
    list_P_size <- parLapply(cl, 3:ncol(df), function(x) { size_P(df[,x]) })
    stopCluster(cl)
    #Used to create empty files to store output lists of breakpoint and cell counts
    #dataframes <- list.files(pattern = "*.txt$")
    #for (i in 1:length(dataframes)) {
    #outputFile <- sub(".txt", "_cell_counts.txt", dataframes[i])
    #file.create(outputFile)}
    
    nb_of_breakpoints[[file]] <- lengths(list_breakpoint_size)
    extremeNb_bp_cells <- which(nb_of_breakpoints[[file]] > quantile(nb_of_breakpoints[[file]], probs = q.lim))
  
    list_breakpoint_coordinates <- list_breakpoint_coordinates[-c(extremeNb_bp_cells)]
    list_breakpoint_counts <- lapply(1:length(list_breakpoint_coordinates), function(x) { breakpoint_counts(list_breakpoint_coordinates[[x]]) }) #changed to remove cells with too many bp
    lapply(1:54, function(x) { breakpoint_counts(list_breakpoint_coordinates[[x]]) }) #changed to remove cells with too many bp
    #Sum list_breakpoint_counts
    bp <- Reduce('+',list_breakpoint_counts[which(lengths(list_breakpoint_counts) > 1)])
    
    #Export the list_breakpoint_counts to corresponding file
    write.table(bp, file = paste0("resolution.95/",sub(".txt", "_bp_counts.txt", dataframes[file])), row.names = FALSE, col.names = FALSE)
    
    
    list_unknown_coordinates <- list_unknown_coordinates[-c(extremeNb_bp_cells)] #added to remove cells with too many bp
    ####---------------
    ###-------------
    #All functions listed below
    ##---------
    #-----
    
    
    #cell_counts <- function(cell) {
    #cell - vector with the unknown coordinates within a cell; 
    #it takes the list_unknown_cooridnate as global variable
    #FUnction only applies to cells with breakpoints
    
    
    cell_counts <- vector(mode = "numeric", length = length(bin_start))
    for (i in 1:length(cell_counts)) {
      cell_counts[i] <- length(list_unknown_coordinates)}
    
    for (cell in 1:length(list_unknown_coordinates)) {
      
      if (length(list_unknown_coordinates[[cell]]) > 1) {
        for (i in seq(1, length(list_unknown_coordinates[[cell]])-1, 2)) {
          
          #Determine bin of the first coordinate
          x <- (max(which(bin_start < list_unknown_coordinates[[cell]][i])))
          #Determine bin of the last coordinate
          y <- (min(which(bin_end > list_unknown_coordinates[[cell]][i+1])))
          
          #If unknown blocks spans more than one bin 
          #Remove one cell from all bins between (non-inclusive)  first and last coordinates of unknown block
          if (x != y) {
            cell_counts[(x+1):(y-1)] <- cell_counts[(x+1):(y-1)] - 1
            #If unknown block is in one bin, counts remain unchanged because will be modified later
          } else { cell_counts <- cell_counts
          }
          
          #Check to see if there is a breakpoint within that cell
          if (length(list_breakpoint_coordinates[[cell]]) > 1) {
            for (coordinate in list_breakpoint_coordinates[[cell]]) {
              #If there is a breakpoint cooridinate that is smaller than first coordinate of unknown and falls in that bin, count is unchanged
              if (coordinate < list_unknown_coordinates[[cell]][i] & coordinate > x) {
                cell_counts[x] <- cell_counts[x]
                #If there is no breakpoint coordinate in that bin, remove one cell from count
              } else {
                cell_counts[x] <- cell_counts[x] - 1
                #If there is a breakpoint coordinate that is bigger than last coordinate of unknown and falls in that bin, count is unchanged
              }
              if (coordinate > list_unknown_coordinates[[cell]][i+1] & coordinate < y ) {
                cell_counts[y] <- cell_counts[y]
                #If there is no breakpoint coordinate in that bin, remove one cell from count
              } else {
                cell_counts[y] <- cell_counts[y] - 1
              }
            }
            #If there are no breakpoints in that cell, remove cell from count
          } else {
            cell_counts[x] <- cell_counts[x] - 1
            cell_counts[y] <- cell_counts[y] - 1
          }
        } 
      }
    }
    
    
    #Export the cell_counts to corresponding file
    cc <- cell_counts
    write.table(cc, file = paste0("resolution.95/", sub(".txt", "_cell_counts.txt", dataframes[file])), row.names = FALSE, col.names = FALSE)
    
    
  } 

####-----------------
##---------------



