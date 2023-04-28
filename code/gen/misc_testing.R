virs <- seq(0, 1000, 1)
c1 <- 50
c2 <- 0.1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 50
c2 <- 1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 0.51
c2 <- 0.1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 0.1
c2 <- 1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')


c1 <- 1
c2 <- 0.6
c3 <- 1e-04
m_m = 1 / 3
R0s <- c()
for (vir in virs) {
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  psi <- 1 / (1 / ((gamma + nat_mort + m_m + mort)) + 4)
  F_matrix = matrix(c(0, beta*N_m/N_m^p, epsilon*beta*N_m/N_m^p,
                      0, 0, 0,
                      0, 0, 0),nrow=3, ncol=3, byrow=T)
  V_matrix = matrix(c(-sigma-nat_mort-m_m, 0, 0,
                      sigma, -m_m-gamma-mort-nat_mort, 0,
                      sigma, 0, -psi),nrow=3, ncol=3, byrow=T)
  V_matrix_inverse = solve(V_matrix)
  G_matrix = F_matrix %*% V_matrix_inverse
  R0_one = abs(eigen(G_matrix)$values[1])
  R0_num = beta * N_m * psi * sigma + epsilon * beta * N_m * (m_m * sigma + gamma * sigma + mort * sigma + sigma * nat_mort)
  R0_denom = N_m^p * (- sigma - nat_mort - m_m) * (- psi) * (- m_m - gamma - mort - nat_mort)
  R0_two = abs(R0_num / R0_denom)
  if (round(R0_one, 2) != round(R0_two, 2)) {
    stop(paste0(R0_one, '; ', R0_two))
  }
  R0s <- c(R0s, R0_one)
}

# Plot
plot(virs, R0s, type='l')
virs[which(R0s == max(R0s))]
#plot(virs * c3, R0s, type='l')



# To include in the testing doc:
# Verified that the carrying capacity now makes sense at DFE when we also exclude arrival of recovereds into market
# To check: is it true that DFE means there are no E, I, and R migrations into markets as well? Since there is no disease in the system? I think so, because there shouldn't be any of these classes at DFE, just S
# To check: is it true that E_f, I_f, and R_f are not disease classes? I think so. If they are not disease classes, and just ways of getting new infections, they will be differentiated away since they do not depend on any of the market compartments. Even if it is, it would only be the exponential assumption that makes it treated as a new incubation, but conceptually, it is not a new infection
# To check: is it true that E_f, I_f, and R_f migrations will not be included in the new infections part of the Greek F matrix? I think so, again, will be differentiated away anyway for some reason as above
