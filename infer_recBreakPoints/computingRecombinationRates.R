#.libPaths( c("/home/isabelalves", .libPaths() ) )
args <- commandArgs(trailingOnly = TRUE)

library(randomcoloR)
chrNb <- as.numeric(args[1])
#chrNb <- 14

print(paste0("Computing recombination rates for chrm: ", chrNb))
chrTag <- paste0("chr", chrNb)
hapMapRecombFile <- paste0("genetic_map_GRCh37_", chrTag, ".txt")
deCodeRecombFile <- paste0("male.map.", chrTag, ".txt")

#path to the centromere file
#/home/ialves/Dropbox/singleCellProject/phasing_donors/preprocessing
openCentLocation <- read.table("/home/ialves/Dropbox/singleCellProject/phasing_donors/preprocessing/centromereOnlyReg.txt", header = F,) 
pathToDepthFiles <- "/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/bulk_depth"

folderPath <- paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/resolution/resolution.95")
lapply(list.files(path=paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag), 
                  pattern = "InferredHapsHMap*"), function(x) { file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/",x),
                                                                          folderPath, recursive = F, copy.mode = T) }) 
file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/", hapMapRecombFile),
          folderPath, recursive = F, copy.mode = T)
file.copy(paste0("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/all_chromosomes_Dec2018/Dec17/", chrTag, "/", deCodeRecombFile),
          folderPath, recursive = F, copy.mode = T)
setwd(folderPath)
require(RColorBrewer)

# [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
# [9] "#999999"
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

#computing recombination rates from deCode
deCode_chr <- data.frame(read.table(file=paste0("male.map.chr", chrNb, ".txt"), header = F, na.strings = "NA"))
deCodecM <- c()

for (binbNb in 1:length(bin_counts)) {
  
  #binNb <- 1
  index_tmp <- which(deCode_chr[,3] >= bin_start[binbNb] & deCode_chr[,3] < bin_end[binbNb])
  
  if (length(index_tmp) > 0) {
    if (sum(is.na(deCode_chr[index_tmp, 4])) != length(deCode_chr[index_tmp, 4])) {
      deCodecM[binbNb] <- sum(deCode_chr[index_tmp, 4][!is.na(deCode_chr[index_tmp, 4])])
    } else {
      deCodecM[binbNb] <- 0
    }
  } else {
    deCodecM[binbNb] <- 0
  }
}
##-----------

#computing HapMap recombination rates
hapMap_chr <- data.frame(read.table(file=paste0("genetic_map_GRCh37_chr", chrNb, ".txt"), header = T, na.strings = "NA"))
hapMapcM <- c()
startRecomb <- 0

for (binNb in 1:length(bin_counts)) {
  
  #binNb <- 1
  index_tmp <- which(hapMap_chr[,2] >= bin_start[binNb] & hapMap_chr[,2] < bin_end[binNb])
  if (length(index_tmp) > 0) {
    
    if (binNb == 1) {
      hapMapcM[binNb] <- max(hapMap_chr[index_tmp,4])  
    } else {
      hapMapcM[binNb] <- max(hapMap_chr[index_tmp,4])-startRecomb
    }
    
    startRecomb <- max(hapMap_chr[index_tmp,4])
  } else {
    
    hapMapcM[binNb] <- 0
  }
  
}

# computing recombination rates from single cell sperm cells
empiricalRecombRates <- list()
QC_tag <- list()
filesToAnalyse <- list.files(pattern = "*bp_counts.txt$")

for (files in 1:length(filesToAnalyse)) {
  #files <- 1
  QC_tag[[files]] <- strsplit(filesToAnalyse[files], split = "_")[[1]][2]
  empiricalCountsF <- scan(file=filesToAnalyse[files], what = numeric())
  empiricalCellCountsF <- scan(file=sub("bp", "cell", filesToAnalyse[files]), what = numeric())
  empiricalRecombRates[[files]] <- (empiricalCountsF/empiricalCellCountsF)*100
  
}
##------------
colors <- distinctColorPalette(length(empiricalRecombRates))

# #Plotting dist of breakpoints per cell 
# for (ind in 1:length(empiricalRecombRates)) {
  
#   if (ind == 1) {
    
#     print(QC_tag[[ind]])
#     plot(density(empiricalRecombRates[[ind]]), col=colors[ind], ylim=c(0,0.5), xlim=c(0,12))
#     print(quantile(empiricalRecombRates[[ind]], probs = 0.95))
#   } else {
#     print(QC_tag[[ind]])
#     lines(density(empiricalRecombRates[[ind]]), col=colors[ind])
#     print(quantile(empiricalRecombRates[[ind]], probs = 0.95))
    
#   }
# }
# empiricalRecombRates[[ind]][which(empiricalRecombRates[[ind]] < quantile(empiricalRecombRates[[ind]], probs = 0.95))]

#########################
##
## computing the cumulative proportions
##
#########################
pdf(file=paste0("genetic_vs_physicalDist_chrm", chrNb, ".pdf"), width=6, height=6)
plot(1:length(hapMapcM), y=cumsum(hapMapcM), type="l", xlim=c(0,max(cumsum(hapMapcM))), ylim=c(0,max(max(cumsum(hapMapcM)), max(unlist(lapply(empiricalRecombRates, cumsum))))), xlab="Physical Distance (Mb)", 
    ylab="Genetic Distance (cM)", col="grey80", lwd=2, lty=4)
lines(1:length(hapMapcM), y=cumsum(deCodecM), type="l", xlim=c(0,max(cumsum(deCodecM))), ylim=c(0,max(cumsum(deCodecM))), xlab="Physical Distance (Mb)", 
ylab="Genetic Distance (cM)", col="black",lwd=2, lty=3)
rect(openCentLocation[which(openCentLocation[,2] == chrTag),3]/1e6, 0, openCentLocation[which(openCentLocation[,2] == chrTag),4]/1e6, max(max(cumsum(hapMapcM)), max(unlist(lapply(empiricalRecombRates, cumsum)))),
        col = "grey88", border = NA)
        
for (files in 1:length(empiricalRecombRates)) {
  
  lines(1:length(hapMapcM), y=cumsum(empiricalRecombRates[[files]]), type="l", xlim=c(0,max(cumsum(deCodecM))), ylim=c(0,max(cumsum(deCodecM))), xlab="Physical Distance (Mb)", 
  ylab="Genetic Distance (cM)", col=colors[files],lwd=2)
   
}
legend("topright", legend = c("deCODE", "HapMap", unlist(QC_tag)), lty = c(3,4,rep(1,length(empiricalRecombRates))), lwd = 2, col = c("black", "grey80", colors), bty = "n", cex = 0.60)
legend("topleft", legend = paste0("chrm", chrNb), bty = "n", cex = 1.5)
dev.off()
##------------
#--------
pdf(file = paste0("recombRates_SNPDensity_chr", chrNb, "_2.pdf"), width = 10, height = 8)
par(mfrow=c(3,1))
##----------
# set up the plot 
y_range <- range(unlist(empiricalRecombRates))
x_range <- range(1:length(bin_counts))
tickMarks <- seq(1,length(bin_counts), by = 20)

plot(x_range, y_range, type="n", xlab=paste0("Chromosome position (Chrm" , chrNb,")"),
     ylab="cM/Mb", xaxt='n')
axis(1, at = tickMarks, labels = format(apply(rbind(bin_start, bin_end), 2, mean), scientific = T, digits = 2)[tickMarks])
legend("topleft", legend = c("deCODE", "HapMap", unlist(QC_tag)), lty = c(3,4,rep(1,length(empiricalRecombRates))), lwd = 2, col = c("black", "grey80", colors), bty = "n", cex = 0.60)
rect(openCentLocation[which(openCentLocation[,2] == chrTag),3]/1e6, 0, openCentLocation[which(openCentLocation[,2] == chrTag),4]/1e6, y_range[2],
     col = "grey88", border = NA)

for (files in 1:length(empiricalRecombRates)) {
  
  lines(unlist(empiricalRecombRates[files]), type = "l", col = colors[files], lwd=1.5)  
    
  
}
lines(deCodecM, type = "l", lty= 3, lwd = 2, col = "black")
lines(hapMapcM, type = "l", lty= 4, lwd = 2, col = "grey80")

##-------------

#Computing SNP density onto another plot
#m10ToAnalyse <- list.files(pattern = "*NbLinks10.txt$")
m5ToAnalyse <- list.files(pattern = paste0("*", chrTag, "_QCvar25_NbCells10_NbLinks5.txt$"))

#m7ToAnalyse <- list.files(pattern = "*NbLinks7.txt$")
mToAnalyse <- sort(m5ToAnalyse)
density_SNP_phase <- list()
mean_SNP_depth <- list()

for (files in 1:length(mToAnalyse)) {

  #files <- 1
  fileName_df <- data.frame(read.table(file = mToAnalyse[files], header = T))
  indID <- unlist(strsplit(mToAnalyse[files], split = "_"))[2]
  openBulkDepth <- read.table(paste0(pathToDepthFiles, "/wgs.", indID, ".", chrTag, ".commonhetSNPs.DP.ldepth"), header = F) 
  
  pos <- fileName_df[,2]
  intersectingPos <- intersect(pos, openBulkDepth[,2])
  print(length(match(intersectingPos, openBulkDepth[,2])))
  sub_bulkDepth <- openBulkDepth[match(intersectingPos, openBulkDepth[,2]),]
  rm(openBulkDepth)
  rm(fileName_df)
  
  dens_tmp <- c()
  dp_tmp <- c()
  for (binNb in 1:length(bin_counts)) {
    #binNb <- 2 
    idx_tmp <- which(pos >= bin_start[binNb] & pos < bin_end[binNb])
    idx_dp_tmp <- which(sub_bulkDepth[,2] >= bin_start[binNb] & sub_bulkDepth[,2] < bin_end[binNb])
    
    if (length(idx_tmp) > 0) { 
      dens_tmp[binNb] <- length(idx_tmp)
      dp_tmp[binNb] <- mean(sub_bulkDepth[idx_dp_tmp,3], na.rm = T)
    } else {
      dens_tmp[binNb] <- 0
      dp_tmp[binNb] <- 0
      
    }
  }
  density_SNP_phase[[files]] <- dens_tmp
  mean_SNP_depth[[files]] <- dp_tmp
}

#plotting SNP density 

#density_SNP_phase[[4]] <- NULL 
# set up the plot 
y_dens_range <- range(unlist(density_SNP_phase))
x_dens_range <- range(1:length(bin_counts))
tickMarks <- seq(1,length(bin_counts), by = 20)


plot(x_dens_range, y_dens_range, type="n", xlab=paste0("Chromosome position (Chrm" , chrNb,")"),
     ylab="SNP number", xaxt='n')
axis(1, at = tickMarks, labels = format(apply(rbind(bin_start, bin_end), 2, mean), scientific = T, digits = 2)[tickMarks])
rect(openCentLocation[which(openCentLocation[,2] == chrTag),3]/1e6, 0, openCentLocation[which(openCentLocation[,2] == chrTag),4]/1e6, y_dens_range[2],
     col = "grey88", border = NA)

#colors <- brewer.pal(length(empiricalRecombRates),"Set1")
#legend("topleft", legend = c("deCODE", "HapMap"), lty = c(3,4), lwd = 2, col = c("black", "grey80"), bty = "n")

for (files in 1:length(density_SNP_phase)) {
  
  lines(unlist(density_SNP_phase[files]), type = "l", col = colors[files], lwd=1.5)  
  
  
}

## plotting mean dp 
y_dp_range <- range(unlist(mean_SNP_depth))
x_dp_range <- range(1:length(bin_counts))
tickMarks <- seq(1,length(bin_counts), by = 20)

plot(x_dp_range, y_dp_range, type="n", xlab=paste0("Chromosome position (Chrm" , chrNb,")"),
     ylab="Mean depth", xaxt='n')
axis(1, at = tickMarks, labels = format(apply(rbind(bin_start, bin_end), 2, mean), scientific = T, digits = 2)[tickMarks])
rect(openCentLocation[which(openCentLocation[,2] == chrTag),3]/1e6, 0, openCentLocation[which(openCentLocation[,2] == chrTag),4]/1e6, y_dp_range[2],
     col = "grey88", border = NA)

for (files in 1:length(mean_SNP_depth)) {
  
  lines(unlist(mean_SNP_depth[files]), type = "l", col = colors[files], lwd=1.5)  
  
  
}


dev.off()
