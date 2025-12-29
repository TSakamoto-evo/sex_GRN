folder <- "max_div/"

ng <- 8

list_sd <- rep(0, ng + ng^2)

for(run in 1:1000){
  if(file.exists(paste0(folder, "folder", run,"/regi_time.txt"))){
    data <- read.table(paste0(folder, "folder", run,"/regi_time.txt"))
    data2 <- read.table(paste0(folder, "folder", run, "/log_association.txt.gz"))
    
    subdata2 <- data2[data2$V1 >= data[1, 1], ]
    subdata2 <- subdata2[order(subdata2$V1)[1], ]
    
    sd <- apply(subdata2[, -1], 1, function(y){
      return(order(y)[1])
    })
    
    for(i in unique(sd)){
      list_sd[i] <- list_sd[i] + sum(sd == i)
    }
  
    rm(data)
    rm(data2)
  }
  
  if(run %% 100 == 0){
    print(run)
  }
}

prop <- list_sd / sum(list_sd)
max_prop <- max(prop)
max_prop <- 0.15

color_func <- colorRamp(c("white", "#0080ff","#00ff80","#ff8000", "#f08080"))

#postscript("tmp.eps", height=3.5, width=3.5)
plot(NA, xlim=c(0, ng + 1.5), ylim=c(0, ng), axes=F, xlab="", ylab="")
for(i in 1:ng){
  index <- i
  col <- rgb(color_func(prop[index] / max_prop) / 256)
  
  polygon(c(ng+0.5, ng+0.5, ng+1.5, ng+1.5), c(ng-i, ng-i+1, ng-i+1, ng-i), 
          col=col, border="black", lwd=1) 
}

for(i in 1:ng){
  for(j in 1:ng){
    index <- ng + (i-1)*ng + j
    col <- rgb(color_func(prop[index] / max_prop) / 256)

    polygon(c(j-1, j-1, j, j), c(ng-i, ng-i+1, ng-i+1, ng-i), 
            col=col, border="black", lwd=1)
  }
}
#dev.off()

#postscript("tmp.eps", height=3.5, width=3.5)
#plot(NA, xlim=c(0, 1), ylim=c(0, 1))
#for(i in 0:100){
#  col <- rgb(color_func(i / 100) / 256)
#  polygon(c(0, 0, 1, 1), c((i-0.5)/100, (i+0.5)/100, (i+0.5)/100, (i-0.5)/100), 
#          col=col, border=NA)
#}
#dev.off()

print(max(prop))
