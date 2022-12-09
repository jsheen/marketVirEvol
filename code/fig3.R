# 1) set function to get R0 in the model ---------------------------------------
get_R0 <- function(beta, N_m, p, epsilon, sigma, nat_mort, m_m, gamma, mort, psi, both=T) {
  R0_num = beta * N_m * psi * sigma + epsilon * beta * N_m * (m_m * sigma + gamma * sigma + mort * sigma + sigma * nat_mort)
  R0_denom = N_m^p * (- sigma - nat_mort - m_m) * (- psi) * (- m_m - gamma - mort - nat_mort)
  R0 = abs(R0_num / R0_denom)
  if (both) {
    F_matrix = matrix(c(0, beta*N_m/N_m^p, epsilon*beta*N_m/N_m^p,
                        0, 0, 0,
                        0, 0, 0),nrow=3, ncol=3, byrow=T)
    V_matrix = matrix(c(-sigma-nat_mort-m_m, 0, 0,
                        sigma, -m_m-gamma-mort-nat_mort, 0,
                        sigma, 0, -psi),nrow=3, ncol=3, byrow=T)
    V_matrix_inverse = solve(V_matrix)
    G_matrix = F_matrix %*% V_matrix_inverse
    R0_mat = abs(eigen(G_matrix)$values[1])
    if (round(R0, 2) != round(R0_mat, 2)) {
      stop(paste0(R0, '; ', R0_mat))
    }
  }
  return(R0)
}

# 2) parameters needed to get R0 for range of psi_cleans and m_ms --------------
c1 = 1
c2 = 0.5
c3 = 1 / 1000
virs = seq(1, 1000, 1)
N_m = 1000
p = 1
epsilon = 0.1
sigma = 1 / 5
nat_mort = 1 / 365
gamma = 1 / 5
psi_cleans <- seq(1, 20, 1)
m_ms <- seq(1/365, 1/3.5, 0.01)
final_ls <- list()
final_ls_dex <- 1
for (psi_clean in psi_cleans) {
  for (m_m in m_ms) {
    R0s <- c()
    for (vir in virs) {
      beta = c1 * (vir) ^ c2
      mort = c3 * vir
      psi = (1 / ((1 / (gamma + nat_mort + m_m + mort)) + 4)) * psi_clean
      R0 <- get_R0(beta=beta, N_m=N_m, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, both=F)
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
ggplot(final_df, aes(psi_clean, m_m)) + geom_tile(aes(fill = opt_vir))


















