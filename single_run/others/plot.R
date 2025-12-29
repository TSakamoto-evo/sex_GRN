data <- read.table("../c++/fitness.txt")

num_g <- 8

postscript("tmp.eps", height=2.8, width=4)
plot(NA, xlim=c(0, 2000000), ylim=c(0, 1), 
     xlab="generation", ylab="proportion")
for(i in 1:(length(data$V1) - 1)){
  t1 <- data$V1[i]
  t2 <- data$V1[i+1]
  
  x1 <- data$V4[i]
  x2 <- data$V4[i+1]
  
  y1 <- x1 + data$V3[i]
  y2 <- x2 + data$V3[i+1]
  
  z1 <- y1 + data$V5[i]
  z2 <- y2 + data$V5[i+1]
  
  polygon(c(t1, t1, t2, t2), c(0, x1, x2, 0), col=2, border=NA)
  polygon(c(t1, t1, t2, t2), c(x1, y1, y2, x2), col=4, border=NA)
  polygon(c(t1, t1, t2, t2), c(y1, z1, z2, y2), col=gray(0.8), border=NA)
}

#abline(v=0)
#abline(v=10000)
#abline(v=50000)
#abline(v=500000)
#abline(v=2000000)

dev.off()
