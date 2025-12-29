folder <- "default/"

#postscript("tmp.eps", height=4, width=4)
#dev.off()

#for(run in 1:1000){
  if(file.exists(paste0(folder, "folder", run,"/log_pval_list.txt"))){
    data <- read.table(paste0(folder, "folder", run,"/log_association.txt.gz"))
    data2 <- read.table(paste0(folder, "folder", run,"/regi_time.txt"))
    data3 <- read.table(paste0(folder, "folder", run,"/log_pval_list.txt"))
    
    data4 <- read.table(paste0(folder, "folder", run,"/fitness.txt.gz"))
    
    max_val <- 200
    
    plot(NA, xlim=c(data2[1, 1]+500000-2000, data2[1, 1]+500000+10000), ylim=c(0, 72), main=run, yaxs="i")
    for(i in 1:72){
      for(j in 1:dim(data)[1]){
        val <- -data[j, i+1] * log10(exp(1)) / max_val
        if(val > 1){
          val <- 1
        }
        col <- rgb(1, 1-val, 1-val)
        
        if(j > 1 & j < dim(data)[1]){
          sep1 <- (data$V1[j] - data$V1[j-1])/2
          sep2 <- (data$V1[j+1] - data$V1[j])/2
        }else if(j==1){
          sep1 <- (data$V1[j+1] - data$V1[j])/2
          sep2 <- sep1
        }else{
          sep1 <- (data$V1[j] - data$V1[j-1])/2
          sep2 <- sep1
        }
        
        polygon(c(data$V1[j]-sep1, data$V1[j]-sep1, data$V1[j]+sep2, data$V1[j]+sep2), 
                c(i, i-1, i-1, i), col=col, border=F)
      }
    }
    abline(v=data2$V1[1]+500000, lwd=1)
    abline(v=data2$V1[2], col=4, lwd=2)
           
    Sys.sleep(3)
  }
  
  #plot(as.vector(as.matrix(-log10(data[, -1]))))
  
#}

#abline(h=56)
