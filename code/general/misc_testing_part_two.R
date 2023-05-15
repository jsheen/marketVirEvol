# Look through the R0s of each orig_extincts_df
orig_extincts_df$R0 <- NA
for (i in 1:nrow(orig_extincts_df)) {
  c1 <- orig_extincts_df$c1[i]
  c2 <- orig_extincts_df$c2[i]
  m_m <- orig_extincts_df$m_m[i]
  psi_clean <- orig_extincts_df$psi_clean[i]
  vir <- orig_extincts_df$opt_vir[i]
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  psi <-  (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
  R0 <- get_R0(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, kappa=psi_clean, both=F)
  orig_extincts_df$R0[i] <- R0
  
  # run sim
  mod_eqn_nob <- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
      dS_m = -(beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m -epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -nat_mort*S_m -m_m*S_m +m_f*S_f
      dE_m = (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m +epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -sigma*E_m -nat_mort*E_m -m_m*E_m  #+m_f*E_f
      dI_m = sigma*E_m -m_m*I_m -gamma*I_m -mort*I_m -nat_mort*I_m #+m_f*I_f
      dR_m = gamma*I_m -m_m*R_m -nat_mort*R_m #+m_f*R_f
      dH_m = sigma*E_m -psi*kappa*H_m #+m_f*E_f
      return(list(c(dS_m, dE_m, dI_m, dR_m, dH_m)))})}
  parameters <- c(m_m = m_m, gamma = gamma, m_f = m_f, beta = beta,
                  S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = nat_mort,
                  p = p, mort = mort, sigma = sigma, epsilon = epsilon, psi=psi, kappa = psi_clean)
  init <- c(S_m = 999,
            E_m = 0,
            I_m = 1,
            R_m = 0,
            H_m = 0)
  times <- seq(1, 730, 1)
  out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
  out.df <- as.data.frame(out)
  if (all(out.df$I_m == cummin(out.df$I_m))) {
    print(i)
  }
}

# Plot some parameter combination to see if it goes extinct
i <- 2
c1 <- orig_extincts_df$c1[i]
c2 <- orig_extincts_df$c2[i]
m_m <- orig_extincts_df$m_m[i]
psi_clean <- orig_extincts_df$psi_clean[i]
opt_vir <- orig_extincts_df$opt_vir[i]
mort <- opt_vir * c3
psi <-  (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
mod_eqn_nob <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS_m = -(beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m -epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -nat_mort*S_m -m_m*S_m +m_f*S_f
    dE_m = (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m +epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -sigma*E_m -nat_mort*E_m -m_m*E_m  #+m_f*E_f
    dI_m = sigma*E_m -m_m*I_m -gamma*I_m -mort*I_m -nat_mort*I_m #+m_f*I_f
    dR_m = gamma*I_m -m_m*R_m -nat_mort*R_m #+m_f*R_f
    dH_m = sigma*E_m -psi*kappa*H_m #+m_f*E_f
    return(list(c(dS_m, dE_m, dI_m, dR_m, dH_m)))})}
beta <- c1 * (opt_vir)^c2
get_R0(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, kappa=psi_clean, both=F)
parameters <- c(m_m = m_m, gamma = gamma, m_f = m_f, beta = beta,
                S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = nat_mort,
                p = p, mort = mort, sigma = sigma, epsilon = epsilon, psi=psi, kappa = psi_clean)
init <- c(S_m = 999,
          E_m = 0,
          I_m = 1,
          R_m = 0,
          H_m = 0)
times <- seq(1, 3650, 1)
out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
out.df <- as.data.frame(out)
plot(out.df$time, out.df$I_m, type='l')

get_R02 <- function(beta, m_f, S_f, b, p, epsilon, sigma, nat_mort, m_m, gamma, mort, psi, kappa, both=T) {
  F_matrix = matrix(c(0, -beta*((m_f * S_f) / (m_m + nat_mort))^(1-p), -epsilon*beta*((m_f * S_f) / (m_m + nat_mort))^(1-p),
                      0, 0, 0,
                      0, 0, 0),nrow=3, ncol=3, byrow=T)
  V_matrix = matrix(c((-sigma - nat_mort - m_m), 0, 0,
                      sigma, (- m_m - gamma - mort - nat_mort), 0,
                      sigma, 0, - kappa * psi), nrow=3, ncol=3, byrow=T)
  V_matrix_inverse = solve(V_matrix)
  G_matrix = F_matrix %*% V_matrix_inverse
  R0_mat = eigen(G_matrix)$values[1]
  return(R0_mat)
}
get_R02(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, kappa=psi_clean, both=F)
