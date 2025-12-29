folder <- "../c++/"

ng <- 8
index <- 48
time <- 200000

cols0 = colorRamp(c("white","black"))
cols1 = colorRamp(c("white","#9400d3"))
cols2 = colorRamp(c("white","#32cd32"))
cols3 = colorRamp(c("white","#ff8c00"))

data <- read.table(paste0(folder, "pca_res.txt"))
data2 <- read.table(paste0(folder, "regi_genotype.txt"))
head(data)

subdata <- data[data$V1 == time, ]
subdata2 <- data2[data2$V1 == time, ]

if(index < ng){
  g1 <- subdata2[, 6+index]
  g2 <- subdata2[, 6+index+ng]
  val_mean <- mean(c(g1, g2))
  geno <- (g1 > val_mean) + (g2 > val_mean)
  
  hist(c(g1, g2), breaks=100)
}else{
  g1 <- subdata2[, 6+2*ng+(index-ng)]
  g2 <- subdata2[, 6+2*ng+ng**2+(index-ng)]
  val_mean <- mean(c(g1, g2))
  geno <- (g1 > val_mean) + (g2 > val_mean)
  
  hist(c(g1, g2), breaks=100)
}

#postscript("tmp.eps", height=3.15, width=2.75)
plot(NA, xlim=c(-max(data$V4), -min(data$V4)), ylim=c(min(data$V3), max(data$V3)), 
     xlab="PC2", ylab="PC1")

geno2 <- numeric(0)
for(i in 1:length(geno)){
  geno2 <- c(geno2, rep(geno[i], max(subdata$V2) + 1))
}

prop <- 0.9 * subdata$V2 / max(subdata$V2) + 0.1
cols <- rgb((cols1(prop) * (geno2 == 0) + cols2(prop) * (geno2 == 1) + cols3(prop) * (geno2 == 2))/256) 

#cols <- rgb(cols0(prop) / 256)

#points(-subdata$V4, subdata$V3, pch=16, cex=1, col=cols)
points(-subdata$V4[geno2 == 2], subdata$V3[geno2 == 2], pch=16, cex=0.33, col=cols[geno2 == 2])
points(-subdata$V4[geno2 == 1], subdata$V3[geno2 == 1], pch=16, cex=0.33, col=cols[geno2 == 1])
points(-subdata$V4[geno2 == 0], subdata$V3[geno2 == 0], pch=16, cex=0.33, col=cols[geno2 == 0])
#dev.off()


#plot(NA, xlim=c(0, 1), ylim=c(0, 1))
#for(i in 0:max(subdata$V2)){
#  prop <- 0.9 * i / max(subdata$V2) + 0.1
#  col1 <- rgb(cols1(prop) / 256)
#  col2 <- rgb(cols2(prop) / 256)
#  col3 <- rgb(cols3(prop) / 256)
#  
#  y_min <- 1-i/max(subdata$V2) - 0.5/max(subdata$V2)
#  y_max <- 1-i/max(subdata$V2) + 0.5/max(subdata$V2)
#  
#  polygon(c(0, 0, 1/3, 1/3), c(y_min, y_max, y_max, y_min), col=col1, border=F)
#  polygon(c(1/3, 1/3, 2/3, 2/3), c(y_min, y_max, y_max, y_min), col=col2, border=F)
#  polygon(c(2/3, 2/3, 1, 1), c(y_min, y_max, y_max, y_min), col=col3, border=F)
#}