# 1) Set function to get R0 in the model when infectious poultry shed constantly to the environment -----
get_R0_shed <- function(beta, m_f, S_f, b, p, epsilon, lambda, sigma, nat_mort, m, gamma, mort, kappa, psi, both=T) {
  R0 = ((m_f * S_f + b) / (m + nat_mort))^(1-p) * (sigma / (sigma + nat_mort + m)) * (beta / (m + gamma + mort + nat_mort) + ((epsilon * beta * lambda) / (kappa * psi * (m + gamma + mort + nat_mort))))
  R0 = abs(R0)
  if (both) {
    F_matrix = matrix(c(0, beta*((m_f * S_f + b) / (m + nat_mort))^(1-p), (epsilon * beta)*((m_f * S_f + b) / (m + nat_mort))^(1-p),
                        0, 0, 0,
                        0, 0, 0),nrow=3, ncol=3, byrow=T)
    V_matrix = matrix(c((sigma + nat_mort + m), 0, 0,
                        -sigma, (m + gamma + mort +nat_mort), 0,
                        0, -lambda, kappa * psi), nrow=3, ncol=3, byrow=T)
    V_matrix_inverse = solve(V_matrix)
    G_matrix = F_matrix %*% V_matrix_inverse
    R0_mat = eigen(G_matrix)$values[1]
    if (round(R0, 2) != round(R0_mat, 2)) {
      stop(paste0(R0, '; ', R0_mat))
    }
  }
  return(R0)
}

# 2) Set function to get R0 in the model when infectious poultry move at different rate -----
get_R0_shed_diff_migrate <- function(beta, m_f, S_f, b, p, epsilon, lambda, sigma, nat_mort, m, gamma, mort, kappa, psi, m_I, both=T) {
  R0 = ((m_f * S_f + b) / (m + nat_mort))^(1-p) * (sigma / (sigma + nat_mort + m)) * (beta / (m_I + gamma + mort + nat_mort) + ((epsilon * beta * lambda) / (kappa * psi * (m_I + gamma + mort + nat_mort))))
  R0 = abs(R0)
  if (both) {
    F_matrix = matrix(c(0, beta*((m_f * S_f + b) / (m + nat_mort))^(1-p), (epsilon * beta)*((m_f * S_f + b) / (m + nat_mort))^(1-p),
                        0, 0, 0,
                        0, 0, 0),nrow=3, ncol=3, byrow=T)
    V_matrix = matrix(c((sigma + nat_mort + m), 0, 0,
                        -sigma, (m_I + gamma + mort +nat_mort), 0,
                        0, -lambda, kappa * psi), nrow=3, ncol=3, byrow=T)
    V_matrix_inverse = solve(V_matrix)
    G_matrix = F_matrix %*% V_matrix_inverse
    R0_mat = eigen(G_matrix)$values[1]
    if (round(R0, 2) != round(R0_mat, 2)) {
      stop(paste0(R0, '; ', R0_mat))
    }
  }
  return(R0)
}




