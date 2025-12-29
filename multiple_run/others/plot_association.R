folder <- "default/"

#postscript("tmp.eps", height=4, width=4)
#dev.off()

run <- 175

#for(run in 1:1000){
  if(file.exists(paste0(folder, "folder", run,"/log_association.txt"))){
    data <- read.table(paste0(folder, "folder", run,"/log_association.txt"))
  }else{
    data <- read.table(paste0(folder, "folder", run,"/log_association.txt.gz"))
  }
  data <- data[data$V1 %% 10000 == 0, ]
  
  max_val <- 200
  
  #plot(NA, xlim=c(2.615e+6, 2.635e+6), ylim=c(0, 72), yaxs="i", main=run)
  plot(NA, xlim=c(0, max(data$V1)), ylim=c(0, 72), yaxs="i", main=run)
  
  for(i in 1:72){
    for(j in 1:dim(data)[1]){
      val <- -data[j, i+1] * log10(exp(1)) / max_val
      if(val > 1){
        val <- 1
      }
      col <- rgb(1, 1-val, 1-val)
      
      if(j == 1){
        xmax <- (data$V1[j+1] - data$V1[j]) / 2
        xmin <- xmax
        
        if(val != 0){
          polygon(c(data$V1[j] - xmin, data$V1[j] - xmin, data$V1[j] + xmax, data$V1[j] + xmax), 
                c(i, i-1, i-1, i), col=col, border=F)
        }
      }else if(j == length(data$V1)){
        xmin <- (data$V1[j] - data$V1[j-1]) / 2
        xmax <- xmin
        
        if(val != 0){
          polygon(c(data$V1[j] - xmin, data$V1[j] - xmin, data$V1[j] + xmax, data$V1[j] + xmax), 
                c(i, i-1, i-1, i), col=col, border=F)
        }
      }else{
        xmin <- (data$V1[j] - data$V1[j-1]) / 2
        xmax <- (data$V1[j+1] - data$V1[j]) / 2
        if(val != 0){
          polygon(c(data$V1[j] - xmin, data$V1[j] - xmin, data$V1[j] + xmax, data$V1[j] + xmax), 
                c(i, i-1, i-1, i), col=col, border=F)
        }
      }
    }
  }
  
  #plot(as.vector(as.matrix(-log10(data[, -1]))))
  
  rm(data)
  Sys.sleep(1)
#}

#abline(h=48)
