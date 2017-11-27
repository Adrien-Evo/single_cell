install.packages("fields")
library("fields")

yCoord <- seq(3, 3*(length(finalInference)+3), by = 3)
pdf(file="segmentDist_Aug22.pdf", height = 12, width = 6)
for (cellNb in 1:50) {
  
  #cellNb <- 2
  if (cellNb == 1) {
    
    plotCell <- rep(0, length(finalInference[[cellNb]]))
    plotCell[which(finalInference[[cellNb]] == "M" | finalInference[[cellNb]] == "P")] <-  yCoord[cellNb]
    colCell[which(finalInference[[cellNb]] == "P")] <- "black"
    colCell[which(finalInference[[cellNb]] == "M")] <- "red"
    colCell[which(finalInference[[cellNb]] == "U" | finalInference[[cellNb]] == ".")] <- "white"
    plot(plotCell, type = "p", pch= "l", col=colCell, cex=0.5, ylim=c(1,(yCoord[(length(yCoord))]+3)))
  } else {
    
    plotCell <- rep(0, length(finalInference[[cellNb]]))
    plotCell[which(finalInference[[cellNb]] == "M" | finalInference[[cellNb]] == "P")] <- yCoord[cellNb]
    colCell[which(finalInference[[cellNb]] == "P")] <- "black"
    colCell[which(finalInference[[cellNb]] == "M")] <- "red"
    colCell[which(finalInference[[cellNb]] == "U" | finalInference[[cellNb]] == ".")] <- "white"
    points(plotCell, type = "p", pch="l", col=colCell, cex=0.5)
    
    
  }
  
  
  
}
dev.off()
# 
# xtmp <- 1:10
# xtmpTwo <- 1:10
# ytmp <-  rep(1,10)
# ytmpTwo <- rep(3,10)
# plot(xtmp, y=ytmp, type = "p", pch= "l", col="red", ylim=c(0,10))
# points(xtmpTwo, y=ytmpTwo, type = "p", pch = "l", col="blue")
