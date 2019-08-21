pathToFile <- "/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17"
Chroms <- paste0("chr", 1:22)

columnsToplot <- c(2:4)
ymaxCol <- c(25, 350, 5000)
ymaxZoomCol <- c(5,170,600)
xlabsCol <- c("Mean nb of bp per cell", "Nb of Cells", "Ave Break point size")
out_fNames <- c("Mean_nb_of_bp_per_cell", "Nb_of_Cells", "Ave_Break_point_size")

for (c in 1:length(columnsToplot)) {

  pdf(file=paste0(pathToFile, "/", out_fNames[c] ,".pdf"), height = 6, width = 9)
  par(mfrow=c(2,1), mar=c(4,4,1,1))
  
  plot(NULL, xlim=c(1,44), ylim=c(0,ymaxCol[c]), ylab=xlabsCol[c], xlab = "Chromosome", xaxt = 'n')
  axis(1, at=seq(1,44,length.out = 22), labels = 1:22)
  xCoor <- seq(1,44,length.out = 22)
    
  COUNT <- 1
  for (chr in Chroms) {
    
    #chr <- "chr1"
    openTableRes <- read.table(file=paste0(pathToFile, "/", chr, "/resolution/resolution.95/summaryBp_", chr, ".txt"), header = T)
    cntrlsCases <- grepl("S" ,as.character(openTableRes[,1]))
    colorsV <- rep("grey50", nrow(openTableRes))
    colorsV[which(cntrlsCases == TRUE)] <- "chocolate1"
    #points(x=rep(1, nrow(openTableRes)), y=openTableRes[,2], pch="", cex=0.5)  
    text(x=rep(xCoor[COUNT], nrow(openTableRes)), y=openTableRes[,columnsToplot[c]], labels = openTableRes[,1], cex=0.75, col = colorsV)
    COUNT <- COUNT+1
  }
  
  plot(NULL, xlim=c(1,44), ylim=c(0,ymaxZoomCol[c]), ylab=xlabsCol[c], xlab = "Chromosome", xaxt = 'n')
  axis(1, at=seq(1,44,length.out = 22), labels = 1:22)
  
  COUNT <- 1
  for (chr in Chroms) {
    
    #chr <- "chr1"
    openTableRes <- read.table(file=paste0(pathToFile, "/", chr, "/resolution/resolution.95/summaryBp_", chr, ".txt"), header = T)
    cntrlsCases <- grepl("S" ,as.character(openTableRes[,1]))
    colorsV <- rep("grey50", nrow(openTableRes))
    colorsV[which(cntrlsCases == TRUE)] <- "chocolate1"
    #points(x=rep(1, nrow(openTableRes)), y=openTableRes[,2], pch="", cex=0.5)  
    text(x=rep(xCoor[COUNT], nrow(openTableRes)), y=openTableRes[,columnsToplot[c]], labels = openTableRes[,1], cex=0.75, col = colorsV)
    COUNT <- COUNT+1
  }
  dev.off()
}
