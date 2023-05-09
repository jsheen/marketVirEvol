# 0) libraries and sources -----------------------------------------------------
library(ggplot2)
source('~/marketVirEvol/code/general/gen_functions.R')

# 1) parameters needed to get R0 for range of psi_cleans and m_ms --------------
c1 = 1
c2 = 0.5
c3 = 1 / 1000
virs = seq(1, 1000, 1)
N_m = 1000
p = 1
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
psi_cleans <- seq(1, 10, 0.05)
m_ms <- seq(1/365, 1/5.5, 0.001)
final_ls <- list()
final_ls_dex <- 1
for (psi_clean in psi_cleans) {
  for (m_m in m_ms) {
    R0s <- c()
    for (vir in virs) {
      beta = c1 * (vir) ^ c2
      mort = c3 * vir
      psi = (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
      R0 <- get_R0(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, kappa=psi_clean, both=F)
      R0s <- c(R0s, R0)
    }
    opt_vir <- virs[which(R0s == max(R0s))]
    new_row <- data.frame(matrix(c(psi_clean, m_m, opt_vir), nrow=1, ncol=3))
    final_ls[[final_ls_dex]] <- new_row
    final_ls_dex <- final_ls_dex + 1
  }
}
final_df <- do.call(rbind, final_ls)
colnames(final_df) <- c('psi_clean', 'm_m', 'opt_vir')

# 3) First, plot the tradeoff curve used ---------------------------------------
betas = c()
morts = c()
for (vir in virs) {
  betas <- c(betas, c1 * vir ^c2)
  morts <- c(morts, c3 * vir)
}
plot(morts, betas, col='red', type='l', lwd=5, main='Transmission-Mortality Tradeoff', xlab='Mortality rate', ylab='Transmission rate')

# 4) Next, plot heatmap --------------------------------------------------------
ggplot(final_df, aes(psi_clean, m_m)) + geom_tile(aes(fill = opt_vir)) + 
  xlab('κ') + ylab(bquote(m[m])) + labs(fill="ESS") + 
  theme(axis.text = element_text(size=20)) + theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_text(size=22)) + theme(legend.text = element_text(size=18))

