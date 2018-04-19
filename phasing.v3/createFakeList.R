nbSamples <- sample(10:100, 10, replace = T)
c_o_exp <- list()
count <- 1

for (nb in nbSamples) {
  
  df_exp <- matrix(sample(0:8, 4*nb, replace = T), ncol = 4, byrow = T)
  
  pos_names_one <- sample(1000:1e5, nrow(df_exp),replace = F)
  pos_names_two <- sample(1000:1e5, nrow(df_exp),replace = F)
  
  df_final <- as.data.frame(cbind(rep("chr1", nrow(df_exp)), pos_names_one, pos_names_two, df_exp))
  names(df_final) <- c("chrNb","PosOne", "PosTwo", "0/0", "0/1", "1/0", "1/1")
  
  
  c_o_exp[[count]] <- subsettingLinksMatrix(df_final, 3)  
  count <- count+1
}
