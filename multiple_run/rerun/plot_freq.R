folder <- "folder73"

data <- read.table(paste0(folder, "/regi_genotype.txt.gz"))
data2 <- read.table(paste0(folder, "/axis_val.txt.gz"))

col1 <- "#ffb6c1"
col2 <- 2
col3 <- 4
col4 <- gray(0.8)

# Focal time range
t_min <- 2326000
t_max <- 2329000

# define threshold genotypic values to classify alleles at each locus
crit1 <- -0.1
crit2 <- 0.05

ng <- 8

# locus indices
index1 <- 40
index2 <- 56

if(index1 < ng){
  # type-G
  list11 <- data[, 6+index1]
  list12 <- data[, 6+ng+index1]
}else{
  # type-B
  list11 <- data[, 6+2*ng+(index1-ng)]
  list12 <- data[, 6+2*ng+ng^2+(index1-ng)]
}

if(index2 < ng){
  # type-G
  list21 <- data[, 6+index2]
  list22 <- data[, 6+ng+index2]
}else{
  # type-B
  list21 <- data[, 6+2*ng+(index2-ng)]
  list22 <- data[, 6+2*ng+ng^2+(index2-ng)]
}

par(mfrow=c(1, 2))
# plot effect-size distribution during the focal time range
# it is useful to determine the threshold value
hist(c(list11, list12), main="locus 1")
hist(c(list21, list22), main="locus 2")
par(mfrow=c(1, 1))

# allele count of each inds
geno1 <- (list11 > crit1) + (list12 > crit1)
geno2 <- (list21 > crit2) + (list22 > crit2)

par(mfrow=c(3, 3))
# check the relationship between genotype and phenotype
for(i in 0:2){
  for(j in 0:2){
    subdata <- data[geno1 == i & geno2 == j, ]
    subdata2 <- data2[geno1 == i & geno2 == j, ]
    
    if(dim(subdata)[1] > 0){
      list <- numeric(0)
      for(k in 1:dim(subdata)[1]){
        list <- c(list, rep(subdata2$V2[k], subdata$V3[k]))
      }
      hist(list, main=paste(i, j))
    }
  }
}
par(mfrow=c(1, 1))

postscript("tmp.eps", height=3.5, width=4)
plot(NA, xlim=c(t_min, t_max), ylim=c(0, 1), xlab="generation", ylab="frequency")

for(i in unique(data$V1)){
  if(i >= t_min & i <= t_max){
    subdata <- data[data$V1 == i, ]
    subgeno1 <- geno1[data$V1 == i]
    subgeno2 <- geno2[data$V1 == i]
    
    # mamually specify the focal genotypes
    num1 <- sum(subdata$V3[subgeno1 == 2 & subgeno2 == 1])
    num2 <- sum(subdata$V3[subgeno1 == 1 & subgeno2 == 0])
    num3 <- sum(subdata$V3[subgeno1 == 2 & subgeno2 == 0])
    
    # other genotypes are colored in gray
    num4 <- sum(subdata$V3) - num1 - num2 - num3
    
    p1 <- num1 / sum(subdata$V3)
    p2 <- (num1 + num2) / sum(subdata$V3)
    p3 <- (num1 + num2 + num3) / sum(subdata$V3)
    
    # frequency is recorded every five generations
    polygon(c(i-5, i-5, i+5, i+5), c(0, p1, p1, 0), col=col1, border=F)
    polygon(c(i-5, i-5, i+5, i+5), c(p1, p2, p2, p1), col=col2, border=F)
    polygon(c(i-5, i-5, i+5, i+5), c(p2, p3, p3, p2), col=col3, border=F)
    polygon(c(i-5, i-5, i+5, i+5), c(p3, 1, 1, p3), col=col4, border=F)
  }
}
dev.off()