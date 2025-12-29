data <- read.table("default/outputs.txt", header=T)
head(data)

list_time2 <- data$t2 - data$t1 - 500001

postscript("tmp.eps", height=2.5, width=4)
hist(log10(list_time2), breaks=20, main="", freq=F)
abline(v=log10(median(list_time2)), col=2, lwd=2)
dev.off()

postscript("tmp2.eps", height=2.5, width=4)
hist(data$W1, breaks=seq(0, 1, 1/50), col="#1B9E77", border=NA, main="", freq=F)
hist(data$W3, breaks=seq(0, 1, 1/50), col="#D95F02", border=NA, add=T, freq=F)
dev.off()

postscript("tmp3.eps", height=2.5, width=4)
hist(data$W2, breaks=seq(0, 1, 1/50), col="#1B9E77", border=NA, main="", freq=F)
hist(data$W4, breaks=seq(0, 1, 1/50), col="#D95F02", border=NA, add=T, freq=F)
dev.off()

list_time3 <- (data$t3 - data$t2)[data$t2 != data$t3]
a <- hist(log10(c(data$t1, list_time3)), breaks=50)

postscript("tmp4.eps", height=2.5, width=4)
hist(log10(list_time3), breaks=a$breaks, col=gray("0.1"), border=NA, freq=F, ylim=c(0, 1.5), main="")
hist(log10(data$t1), breaks=a$breaks, col=gray("0.8"), border=NA, freq=F, add=T)
dev.off()

sum(data$t2 == data$t3)
length(list_time3)

mean(data$t1)

sum(data$old == data$new)
sum(data$old != data$new)

sum(list_time2 <= 50)
