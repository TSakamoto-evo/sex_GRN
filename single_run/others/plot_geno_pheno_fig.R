data <- read.table("../c++/regi_genotype.txt")
head(data)
num_g <- 8

max_g <- 10
max_b <- max(abs(data[,-(1:(5+2*num_g))]))
#max_b <- 0.44

max_geno <- 6
dev_step <- 40

time <- 200000

subdata <- data[data[,1] == time, ]
subdata <- subdata[order(subdata[,3], decreasing=T),]
ng <- dim(subdata)[1]

regi_fitm <- numeric(0)
regi_fitf <- numeric(0)

for(row in 1:min(max_geno, ng)){
  postscript(paste("tmp", row, ".eps", sep=""), width=6.4, height=3.9)
  plot(NA, xlim=c(0, 2*num_g + 5.5), ylim=c(0, num_g+1), axes=F, xlab="", ylab="")
  
  full <- subdata[row,]
  
  regi_fitm <- c(regi_fitm, full[4])
  regi_fitf <- c(regi_fitf, full[5])
  
  full <- full[-1:-5]
  
  g1 <- full[1:num_g]; full <- full[-(1:num_g)]
  g2 <- full[1:num_g]; full <- full[-(1:num_g)]
  b1 <- full[1:(num_g**2)]; full <- full[-(1:(num_g**2))];
  b2 <- full
  
  b1 <- matrix(as.numeric(b1), num_g, num_g, byrow=T)
  b2 <- matrix(as.numeric(b2), num_g, num_g, byrow=T)
  
  p_vec <- as.numeric(g1 + g2)
  for(i in 1:dev_step){
    eff1 <- b1 %*% p_vec
    eff2 <- b2 %*% p_vec
    p_vec <- 0.8 * p_vec + 0.5 + 0.5 * tanh(eff1) + 0.5 + 0.5 * tanh(eff2)
  }
  
  for(i in 1:num_g){
    for(j in 1:num_g){
      if(b1[i, j] > 0){
        colb1 <- rgb(b1[i, j]/max_b, 0, 0)
      }else{
        colb1 <- rgb(0, -b1[i, j]/max_b, 0)
      }
      
      if(b2[i, j] > 0){
        colb2 <- rgb(b2[i, j]/max_b, 0, 0)
      }else{
        colb2 <- rgb(0, -b2[i, j]/max_b, 0)
      }
      
      polygon(c(j, j, j-1, j-1), c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), 
              col=colb1, border="white", lwd=0.5)
      polygon(c(num_g+j+2, num_g+j+2, num_g+j+1, num_g+j+1), 
              c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), col=colb2, border="white", lwd=0.5)
    }
    
    colg1 <- rgb(g1[i]*2/max_g, 0, 0)
    colg2 <- rgb(g2[i]*2/max_g, 0, 0)
    colp <- rgb(p_vec[i]/max_g, 0, 0)
    
    polygon(c(num_g+0.5, num_g+0.5, num_g+1.5, num_g+1.5), 
            c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), col=colg1, border="white", lwd=0.5)
    polygon(c(2*num_g+2.5, 2*num_g+2.5, 2*num_g+3.5, 2*num_g+3.5), 
            c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), col=colg2, border="white", lwd=0.5)
    polygon(c(2*num_g+4.5, 2*num_g+4.5, 2*num_g+5.5, 2*num_g+5.5), 
            c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), col=colp, border="white", lwd=0.5)
  }
  dev.off()
}

par(mfrow=c(1,1))
