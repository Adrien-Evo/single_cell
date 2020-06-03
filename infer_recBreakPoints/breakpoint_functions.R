#Computes breakpoints in every cell of file
breakpoint_computation <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,161]
  
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
  
  #Breakpoint computation requires min. 2 values to compute breakpoint between M and P
  if(length(bp) > 1) {
    
    #Sorts coordinates in ascending order
    bp <- sort(c(bp))
    #print(orderedbp)
    
    #Stores size of breakpoints, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the last and first coordinates of different M/P blocks  
    bpBlocks <- c()
    #If adjacent coordinates in bp are different states (M and P), compute the breakpoint size 
    for (i in 1:(length(bp)-1)) {
      if (v[bp[i+1]] != v[bp[i]]) {
        bpBlocks <- c(bpBlocks, (df[,2][bp[i+1]] - df[,2][bp[i]]))
      }
    }
    #No breakpoint computation if there is only one parental state in cell
  } else {
    
    bpBlocks <- c() #changed Dec 2018
  }
  
  #Returns the size of different breakpoints 
  #print(bpBlocks) 
  return(bpBlocks)
}
##----
#--

#Function to compute size of unknown blocks
size_unknown <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
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
  
  
  #Removes all NA/null values of uCoord
  uCoord <- uCoord[is.finite(uCoord)]
  
  #size_unknown computation requires min. 2 values to compute distance between first and last U coordinate
  if(length(uCoord) > 1) {
    
    #Sorts coordinates of uCoord in ascending order
    #uCoord <- sort(c(uCoord))
    #print(ordereduCoord)
    
    #Stores size of unknown blocks, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the first and last coordinates of U blocks  
    uCoordBlocks <- c()
    #Compute difference between last and first point of U blocks
    for (i in seq(1, (length(uCoord)-1), 2)) {
      uCoordBlocks <- c(uCoordBlocks, (df[,2][uCoord[i+1]] - df[,2][uCoord[i]]))
      
    }
    #No size_unknown computation if there are no U blocks
  } else {
    
    uCoordBlocks <- 0
  }
  
  #Returns the size of different unknown blocks 
  #print(uCoordBlocks) 
  return(uCoordBlocks)
}
##---------
#-----


#Function to compute size of M blocks
size_M <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,109]
  
  #mCoord stores first and last coordinates of M blocks
  mCoord <- c()
  #Neutral parental state of the chromosome
  state <- "I"
  
  #Scans cell from  first to last position of chromosome to detect M positions
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of new M block
    if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores first cooridnate of M block 
      mCoord <- c(mCoord, i) 
      
      #Scans for new U block and finds previous end point of M block
    } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
      state <- "U"
      #Stores M block end coordinate
      mCoord <- c(mCoord, suppressWarnings(max(which(v[1:i] == "M"))))
      
      #Scans for new P block and finds previous end point of M block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores M block end cooridnate in mCoord
      mCoord <- c(mCoord, suppressWarnings(max(which(v[1:i] == "M"))))
      
      #If M is the last state of chromosome, finds end point of M block
    } else if (i == length(v) & state == "M") {
      mCoord <- c(mCoord, suppressWarnings(max(which(v == "M"))))
    }
    
  } #end of going along every element of the vector with the parental states
  
  
  #Removes all NA/null values of mCoord
  mCoord <- mCoord[is.finite(mCoord)]
  
  #size_M computation requires min. 2 values to compute distance between first and last M coordinate
  if(length(mCoord) > 1) {
    
    #Sorts coordinates of mCoord in ascending order
    #mCoord <- sort(c(mCoord))
    #print(mCoord)
    
    #Stores size of M blocks, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the first and last coordinates of M blocks  
    mCoordBlocks <- c()
    #Compute difference between last and first point of M blocks
    for (i in seq(1, (length(mCoord)-1), 2)) {
      mCoordBlocks <- c(mCoordBlocks, (df[,2][mCoord[i+1]] - df[,2][mCoord[i]]))
      
    }
    #No size_unknown computation if there are no M blocks
  } else {
    
    mCoordBlocks <- 0
  }
  
  #Returns the size of different M blocks 
  #print(mCoordBlocks) 
  return(mCoordBlocks)
}
##------------
#--------

#Function to compute size of P blocks
size_P <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,122]
  
  #pCoord stores first and last coordinates of P blocks
  pCoord <- c()
  #Neutral parental state of the chromosome
  state <- "I"
  
  #Scans cell from  first to last position of chromosome to detect P positions
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of new M block
    if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores first cooridnate of P block 
      pCoord <- c(pCoord, i) 
      
      #Scans for new U block and finds previous end point of P block
    } else if (v[i] != "." & v[i] == "U" & v[i] != state) {
      state <- "U"
      #Stores P block end coordinate
      pCoord <- c(pCoord, suppressWarnings(max(which(v[1:i] == "P"))))
      
      #Scans for new M block and finds previous end point of P block
    } else if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores P block end cooridnate in pCoord
      pCoord <- c(pCoord, suppressWarnings(max(which(v[1:i] == "P"))))
      
      #If P is the last state of chromosome, finds end point of P block
    } else if (i == length(v) & state == "P") {
      pCoord <- c(pCoord, suppressWarnings(max(which(v == "P"))))
    }
    
  } #end of going along every element of the vector with the parental states
  
  
  #Removes all NA/null values of pCoord
  pCoord <- pCoord[is.finite(pCoord)]
  
  #size_P computation requires min. 2 values to compute distance between first and last P coordinate
  if(length(pCoord) > 1) {
    
    #Sorts coordinates of pCoord in ascending order
    #pCoord <- sort(c(pCoord))
    #print(pCoord)
    
    #Stores size of P blocks, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the first and last coordinates of P blocks  
    pCoordBlocks <- c()
    #Adjacent coordinates represent first and last point of P block
    for (i in seq(1, (length(pCoord)-1), 2)) {
      pCoordBlocks <- c(pCoordBlocks, (df[,2][pCoord[i+1]] - df[,2][pCoord[i]]))
      
    }
    #No size_unknown computation if there are no P blocks
  } else {
    
    pCoordBlocks <- 0
  }
  
  #Returns the size of different P blocks 
  #print(pCoordBlocks) 
  return(pCoordBlocks)
}

