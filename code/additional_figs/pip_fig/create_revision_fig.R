library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(treemapify)

plot_ls <- list()

# Second, low kappa
pip_toPlot <- read.csv('~/Desktop/stepSize10/kapLowMat2_corrected.csv')
colnames(pip_toPlot) <- seq(1, 1000, 10)
rownames(pip_toPlot) <- rev(seq(1, 1000, 10))
pip_toPlot$id = rev(seq(1, 1000, 10))
melted <- melt(pip_toPlot, id.var='id')
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
melted$value <- ifelse(melted$value == 2, NA, melted$value)
melted$resident_virulence <- as.numeric(as.character(melted$resident_virulence))
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("white", "black"), na.value='black') +theme(legend.position="none", axis.text=element_text(size=15),
                                                                           axis.title=element_text(size=20),
                                                                           plot.title=element_text(size=20,face="bold")) +
  xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('PIP, Low cleaning') + scale_x_continuous(breaks = seq(0, 1000, by = 250))
temp_plot
plot_ls[[1]] <- temp_plot

# Second, high kappa
pip_toPlot <- read.csv('~/Desktop/stepSize10/kapHighMat2_corrected.csv')
colnames(pip_toPlot) <- seq(1, 1000, 10)
rownames(pip_toPlot) <- rev(seq(1, 1000, 10))
pip_toPlot$id = rev(seq(1, 1000, 10))
melted <- melt(pip_toPlot, id.var='id')
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
melted$value <- ifelse(melted$value == 2, NA, melted$value)
melted$resident_virulence <- as.numeric(as.character(melted$resident_virulence))
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("white", "black"), na.value='black') +theme(legend.position="none", axis.text=element_text(size=15),
                                                                           axis.title=element_text(size=20),
                                                                           plot.title=element_text(size=20,face="bold")) +
  xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('PIP, High cleaning') + scale_x_continuous(breaks = seq(0, 1000, by = 250))
temp_plot
plot_ls[[2]] <- temp_plot

# Second, low kappa
pip_toPlot <- read.csv('~/Desktop/stepSize10/kapLowMat2_corrected.csv')
colnames(pip_toPlot) <- seq(1, 1000, 10)
rownames(pip_toPlot) <- rev(seq(1, 1000, 10))
pip_toPlot$id = rev(seq(1, 1000, 10))
melted <- melt(pip_toPlot, id.var='id')
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
melted$value <- ifelse(melted$value == 2, NA, melted$value)
melted$resident_virulence <- as.numeric(as.character(melted$resident_virulence))
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none", axis.text=element_text(size=15),
                                                                          axis.title=element_text(size=20),
                                                                          plot.title=element_text(size=20,face="bold")) +
  xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('PIP + erad. res., Low cleaning') + scale_x_continuous(breaks = seq(0, 1000, by = 250))
temp_plot
plot_ls[[3]] <- temp_plot

# Second, high kappa
pip_toPlot <- read.csv('~/Desktop/stepSize10/kapHighMat2_corrected.csv')
colnames(pip_toPlot) <- seq(1, 1000, 10)
rownames(pip_toPlot) <- rev(seq(1, 1000, 10))
pip_toPlot$id = rev(seq(1, 1000, 10))
melted <- melt(pip_toPlot, id.var='id')
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
melted$value <- ifelse(melted$value == 2, NA, melted$value)
melted$resident_virulence <- as.numeric(as.character(melted$resident_virulence))
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none", axis.text=element_text(size=15),
                                                                          axis.title=element_text(size=20),
                                                                          plot.title=element_text(size=20,face="bold")) +
  xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('PIP + erad. res., High cleaning') + scale_x_continuous(breaks = seq(0, 1000, by = 250))
temp_plot
plot_ls[[4]] <- temp_plot

ggsave(filename=paste0("~/Desktop/out.png"), marrangeGrob(grobs = plot_ls, nrow=2, ncol=2, top=NULL, common.legend = FALSE), width=10, height=10, units='in', dpi=600)


# 
# # First, low kappa
# pip_toPlot <- read.csv('~/Desktop/stepSize100/kapLowMat.csv')
# pip_toPlot <- pip_toPlot[,c(2:11)]
# pip_toPlot[6,8] <- 1
# colnames(pip_toPlot) <- seq(1, 1000, 100)
# rownames(pip_toPlot) <- rev(seq(1, 1000, 100))
# pip_toPlot$id = rev(seq(1, 1000, 100))
# melted <- melt(pip_toPlot, id.var='id')
# colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
# melted$value <- ifelse(melted$value == 2, NA, melted$value)
# temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
#   scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none") +
#   xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('Set 1: low κ')
# temp_plot
# plot_ls[[1]] <- temp_plot
# 
# # First, high kappa
# pip_toPlot <- read.csv('~/Desktop/stepSize100/kapHighMat.csv')
# pip_toPlot <- pip_toPlot[,c(2:11)]
# pip_toPlot[3,7] <- 1
# colnames(pip_toPlot) <- seq(1, 1000, 100)
# rownames(pip_toPlot) <- rev(seq(1, 1000, 100))
# pip_toPlot$id = rev(seq(1, 1000, 100))
# melted <- melt(pip_toPlot, id.var='id')
# colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
# melted$value <- ifelse(melted$value == 2, NA, melted$value)
# temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
#   scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none") +
#   xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('Set 1: high κ')
# temp_plot
# plot_ls[[2]] <- temp_plot
# 
# # Third, low kappa
# pip_toPlot <- read.csv('~/Desktop/stepSize100/kapLowMat3.csv')
# pip_toPlot <- pip_toPlot[,c(2:11)]
# pip_toPlot[4,8] <- 1
# colnames(pip_toPlot) <- seq(1, 1000, 100)
# rownames(pip_toPlot) <- rev(seq(1, 1000, 100))
# pip_toPlot$id = rev(seq(1, 1000, 100))
# melted <- melt(pip_toPlot, id.var='id')
# colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
# melted$value <- ifelse(melted$value == 2, NA, melted$value)
# temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
#   scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none") +
#   xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('Set 3: low κ')
# temp_plot
# plot_ls[[5]] <- temp_plot
# 
# # Third, high kappa
# pip_toPlot <- read.csv('~/Desktop/stepSize100/kapHighMat3.csv')
# pip_toPlot <- pip_toPlot[,c(2:11)]
# pip_toPlot[2,10] <- 1
# colnames(pip_toPlot) <- seq(1, 1000, 100)
# rownames(pip_toPlot) <- rev(seq(1, 1000, 100))
# pip_toPlot$id = rev(seq(1, 1000, 100))
# melted <- melt(pip_toPlot, id.var='id')
# colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
# melted$value <- ifelse(melted$value == 2, NA, melted$value)
# temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
#   scale_fill_manual(values = c("white", "black"), na.value='gray') +theme(legend.position="none") +
#   xlab('Resident virulence') + ylab('Invader virulence') + ggtitle('Set 3: low κ')
# temp_plot
# plot_ls[[6]] <- temp_plot
# 
# 
# 
