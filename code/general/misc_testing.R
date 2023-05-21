# Freq. dependent
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

# Dens. dependent
virs <- seq(0, 1000, 1) # Should infect at least one per day
DFE_markets <- 2267.22
c2 <- 0.1
c1 <- 1 / 45
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect no more than 200 poultry initially
DFE_markets <- 2267.22
c2 <- 0.1
c1 <- 1 / 22
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect at least one per day
DFE_markets <- 2267.22
c2 <- 1
c1 <- 1 / 2300000
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect no more than 200 poultry initially
DFE_markets <- 2267.22
c2 <- 1
c1 <- 1 / 11300
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)



# Else
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

