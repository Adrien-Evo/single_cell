subsettingLinksMatrix <- function(PC_tbl, r) { #local variable: PC_tbl = matrix with HQ link; global variables: 1) minNbLinks; FUNCTIONS: 1) countingHetLinks.
  #r = ratio 
  #PC_tbl <- fullHQlinksComb
  #keep pwise Comb with link counts > minNbLinks
  linksM <- lapply(1:nrow(PC_tbl[,4:7]), function(x) { which(as.numeric(as.matrix(PC_tbl[x,4:7])) >= minNbLinks) })
  vvv <- which(unlist(lapply(linksM, FUN=length)) == 2) #vector containing those positions with two most prevalent links
  
  #write.table(kkk[vvv,], file="clean_AD393.chr15.pairwiseComb.phase2.txt", quote = F, row.names = F, col.names = F)
  
  #subsetting the original table with pwise Comb counts
  matrix_counts_links_ltTWO <- PC_tbl[vvv,]
  
  #computing a table with the link count order 2,3,4,1 means that the hights nb of counts if in column 2 and 3
  order_m <- t(apply(matrix_counts_links_ltTWO[,4:7],1, function(o) { order(o, decreasing = T) })) #order of the most prevalent links among sites in vvv
  #getting the row indexes compatible with SNP het status
  tmpOneColHet <- countingHetLinks(order_m)
  tmpTwoColHet <- which(lapply(1:nrow(matrix_counts_links_ltTWO), function(x) { sum(as.numeric(as.matrix(matrix_counts_links_ltTWO[x,order_m[x,1:2]+3]))) > r*sum(as.numeric(as.matrix(matrix_counts_links_ltTWO[x,order_m[x,3:4]+3]))) }) == T)
  colHet <- intersect(tmpOneColHet,tmpTwoColHet)
  
  #building the m_hap matrix, it includes all the haplotypes that and needs to be clean to rm duplicates and merge overlapping haplotypes
  subMatrix <- matrix_counts_links_ltTWO[colHet,]
  #match(sort(subMatrix[1,3:6], decreasing=T), subMatrix[1,3:6])
  subOrder <- order_m[colHet,]
  
  return(list(subMatrix, subOrder))
  
}
##--------------------
#--------------

#this function takes the HQ filtered matrix with the link counts order and 
#checks whether the two highest links indicate a hetSNP AS THEY SHOULD
countingHetLinks <- function(linkCountM) { #it takes the order matrix and it returns a vector with row indexes
  
  colHet <- c()
  
  for (col in 1:nrow(linkCountM)) {
    
    index_Supported_Links <- linkCountM[col,1:2]
    #print(index_Supported_Links)
    
    if (is.element(1,index_Supported_Links) & is.element(4,index_Supported_Links)) {
      
      #      print("Both first and second pos are heterozygous")
      colHet <- c(colHet,col)
      
    } else if (is.element(2,index_Supported_Links) & is.element(3,index_Supported_Links)) {
      
      #      print("Both first and second pos are heterozygous") 
      colHet <- c(colHet,col)
      
    } 
    
  } 
  return(colHet)
}
##--------------------
#--------------