library(patchwork)
library(ggplot2)

folders <- c("small_mu_g/", "default/", "large_mu_g/", 
             "small_mu_b/", "default/", "large_mu_b/", 
             "large_sigma/", "default/", "small_sigma/",
             "default/", "max_div/")

col1 <- "#cccccc"
col2 <- "#969696"
col3 <- "#525252"
col_list <- c(rep(c(col1, col2, col3), 3), col2, col3)

list_time <- numeric(0)
order <- c()
count <- numeric(0)
no_to <- numeric(0)
count_order <- c()

for(f_index in 1:length(folders)){
  folder <- folders[f_index]
  
  tmp_count <- 0
  for(run in 1:1000){
    if(file.exists(paste0(folder, "folder", run,"/regi_time.txt"))){
      data <- read.table(paste0(folder, "folder", run,"/regi_time.txt"))
      list_time <- c(list_time, log10(data[1, 1]))
      order <- c(order, LETTERS[f_index])
      rm(data) 
    }else{
      tmp_count <- tmp_count + 1
    }
  }
  count <- c(count, tmp_count)
  count_order <- c(count_order, LETTERS[f_index])
  
  data2 <- read.table(paste0(folder, "/turnovers_list.txt"))
  no_to <- c(no_to, sum(data2$V2 > 0))
  
  rm(data2)
}

df <- data.frame(time = list_time, order = order)
dfp <- data.frame(order = count_order, prop = 1-count/1000)
dft <- data.frame(order = count_order, turnover = no_to/(1000 - count))

g1 <- ggplot(dfp, aes(x = order, y = prop, fill=order)) +
  geom_bar(stat="identity", width = 0.5) +
  scale_fill_manual(values = col_list) +
  scale_y_continuous(breaks=c(0, 1), expand = c(0, 0), 
                     limits = c(0, 1)) +
  theme(axis.text.x  = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "none")

g2 <- ggplot(df, aes(x = order, y = time, fill=order))
g2 <- g2 + geom_violin(linewidth = 0.2) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.5,
               colour = "red")+
  scale_fill_manual(values = col_list) +
  scale_y_continuous(limits = c(min(list_time, 4), max(list_time))) + 
  theme(axis.text.x  = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "none")

g3 <- ggplot(dft, aes(x = order, y = turnover, fill=order)) +
  geom_bar(stat="identity", width = 0.5) +
  scale_fill_manual(values = col_list) +
  scale_y_continuous(breaks=c(0, 1), expand = c(0, 0), 
                     limits = c(0, 1)) +
  theme(axis.text.x  = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "none")

postscript("tmp.eps", width=6, height=6)
g3 + g2 + g1 + plot_layout(ncol = 1, height = c(1, 3, 1))
dev.off()
