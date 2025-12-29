library("vioplot")

folders <- c("small_mu_g", "default", "large_mu_g", 
             "small_mu_b", "default", "large_mu_b", 
             "large_sigma", "default", "small_sigma", 
             "default", "max_div")

at_list = c(1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14)

df <- read.table(paste0(folders[1], "/outputs.txt"), header=T)
df$para <- 1

for(i in 2:length(folders)){
  df_add <- read.table(paste0(folders[i], "/outputs.txt"), header=T)
  df_add$para <- i
  df <- rbind(df, df_add)
}


## Time until removal
postscript("tmp.eps", height=3, width=6)
vioplot(log10(df$t2-500001-df$t1)~df$para, ylim=c(0, 7), 
        drawRect = FALSE, col=gray(0.6), at=at_list)

for (i in 1:length(folders)) {
  val <- log10(median((df$t2-500001-df$t1)[df$para == i]))
  segments(x0 = at_list[i] - 0.3, x1 = at_list[i] + 0.3, y0 = val, y1 = val,
           col = "black", lwd = 2)
}
dev.off()

## Heterogametic
postscript("tmp2.eps", height=3, width=6)
vioplot(df$W1~df$para, ylim=c(0, 1), 
        drawRect = FALSE, col="#1B9E77", side = "right", at=at_list, lwd=0.5)

vioplot(df$W3~df$para, 
        drawRect = FALSE, col="#D95F02", side = "left", at=at_list, lwd=0.5, add=T)
dev.off()

## Homogametic
postscript("tmp3.eps", height=3, width=6)
vioplot(df$W2~df$para, ylim=c(0, 1), 
        drawRect = FALSE, col="#1B9E77", side = "right", at=at_list, lwd=0.5)

vioplot(df$W4~df$para, 
        drawRect = FALSE, col="#D95F02", side = "left", at=at_list, lwd=0.5, add=T)
dev.off()

## Time until re-establishment given the fitness reduction
postscript("tmp4.eps", height=3, width=6)
df2 <- df[df$t2 != df$t3, ]
vioplot(log10(df2$t3 - df2$t2)~df2$para, 
        drawRect = FALSE, col=gray("0.1"), side = "left", at=at_list, lwd=0.5, ylim=c(0, 7))
vioplot(log10(df$t1)~df$para, 
        drawRect = FALSE, col=gray("0.8"), side = "right", at=at_list, lwd=0.5, add=T)
usr <- par("usr")
dev.off()

## Turnover probability
postscript("tmp5.eps", height=3, width=6)
prop <- numeric(0)
for(i in 1:length(folders)){
  subdf <- df[df$para == i, ]
  same <- sum(subdf$old == subdf$new)
  diff <- sum(subdf$old != subdf$new)
  
  prop <- c(prop, diff / (same + diff))
  print(same + diff)
}

plot(NA, xlim = usr[1:2], ylim = c(0, 1), xlab = "", ylab = "", xaxt="n", yaxs="i", xaxs="i")
axis(1, at=at_list)
for(i in 1:length(folders)){
  polygon(c(at_list[i] - 0.3, at_list[i] - 0.3, at_list[i] + 0.3, at_list[i] + 0.3), 
          c(0, prop[i], prop[i], 0), col=gray(0.5))
}
dev.off()

## Turnover probability
postscript("tmp6.eps", height=3, width=6)
subdf <- df[df$al1 - df$intro > 0, ]
vioplot(log10(subdf$al1 - subdf$intro)~subdf$para, ylim=c(0, 7), 
        drawRect = FALSE, col=gray("0.8"), side = "right", at=at_list, lwd=0.5)

vioplot(log10(subdf$pgd - subdf$al1)~subdf$para, 
        drawRect = FALSE, col=gray("0.1"), side = "left", at=at_list, lwd=0.5, add=T)
dev.off()


## NA probability
postscript("tmp7.eps", height=3, width=6)
prop <- numeric(0)
for(i in 1:length(folders)){
  subdf <- df[df$para == i, ]
  run_nonna <- sum(is.na(subdf$t3))
  run_na <- sum(!is.na(subdf$t3))
  
  prop <- c(prop, run_nonna / (run_nonna + run_na))
  print(run_nonna)
}

plot(NA, xlim = usr[1:2], ylim = c(0, 1), xlab = "", ylab = "", xaxt="n", yaxs="i", xaxs="i")
axis(1, at=at_list)
for(i in 1:length(folders)){
  polygon(c(at_list[i] - 0.3, at_list[i] - 0.3, at_list[i] + 0.3, at_list[i] + 0.3), 
          c(0, prop[i], prop[i], 0), col=gray(0.5))
}
dev.off()