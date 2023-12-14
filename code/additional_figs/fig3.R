# 0) libraries and sources -----------------------------------------------------
library(ggplot2)
library(latex2exp)
library(RColorBrewer)
source('~/marketVirEvol/code/general/gen_functions.R')

# 1) parameters needed to get R0 for range of psi_cleans and ms --------------
c1 = 1 / 2300
c2 = 0.6
c3 = 1 / 1000
virs = seq(1, 1000, 1)
N = 1000
p = 0
b=0
psi = 1 / 5
m_f = 0.1 / 120
N_f = 1e6
prev_f = 0.12
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.402
R_f = seroprev_f * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f
epsilon = 1
sigma = 1 / 5
nat_mort = 1 / 365
gamma = 1 / 5
phi = 1
psi_cleans <- seq(1, 10, 0.1)
ms <- seq(1/30, 1/5.5, 0.01)
final_ls <- list()
final_ls_dex <- 1
for (psi_clean in psi_cleans) {
  for (m in ms) {
    R0s <- c()
    for (vir in virs) {
      beta = c1 * (vir) ^ c2
      mort = c3 * vir
      lambda = beta * phi
      R0 <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda, 
                        sigma=sigma, nat_mort=nat_mort, m=m, gamma=gamma, mort=mort, kappa=psi_clean, psi=psi, both=F)
      R0s <- c(R0s, R0)
    }
    opt_vir <- virs[which(R0s == max(R0s))]
    new_row <- data.frame(matrix(c(psi_clean, m, opt_vir, max(R0s)), nrow=1, ncol=4))
    final_ls[[final_ls_dex]] <- new_row
    final_ls_dex <- final_ls_dex + 1
  }
}
final_df <- do.call(rbind, final_ls)
colnames(final_df) <- c('psi_clean', 'm', 'opt_vir', 'max_R0')

# 3) First, plot the tradeoff curve used ---------------------------------------
betas = c()
morts = c()
for (vir in virs) {
  betas <- c(betas, c1 * vir ^c2)
  morts <- c(morts, c3 * vir)
}
#plot(morts, betas, col='red', type='l', lwd=5, main='Transmission-Mortality Tradeoff', xlab='Mortality rate', ylab='Transmission rate')

# 4) Next, plot heatmap of optimal virulences ----------------------------------
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(100))
p <- ggplot(final_df, aes(psi_clean, m)) + geom_tile(aes(fill = opt_vir)) + 
  xlab(TeX("$\\kappa$"))  + ylab(expression(italic('m'))) + labs(fill="ESS") +  sc+
  theme(axis.text = element_text(size=8)) + theme(axis.title = element_text(size=10)) +
  theme(legend.title = element_text(size=10)) + theme(legend.text = element_text(size=10)) #+
  #theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10, unit = "pt"))

ggsave("~/Desktop/Fig3.pdf", p, width=3.46, height=2.5, 
        units = "in", device = pdf)

# 5) Next, plot heatmap of optimal R0s (note, this is using the market DFE, and would change under farm conditions) ---------
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(100))
ggplot(final_df, aes(psi_clean, m)) + geom_tile(aes(fill = max_R0)) + 
  xlab(TeX("$\\kappa$"))  + ylab(expression(italic('m'))) + labs(fill="R0") +  sc+
  theme(axis.text = element_text(size=20)) + theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_text(size=22)) + theme(legend.text = element_text(size=18))

