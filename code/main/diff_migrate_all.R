# ------------------------------------------------------------------------------
# @description: this is used to evaluate the first part of the virulence 
#               evolution of markets paper. It is divided into two parts: 
#               - (a) all questions besides extinction question, which requires simulation
#               - (b) extinction question, which requires simulation
# ------------------------------------------------------------------------------
# 0) libraries and sources -----------------------------------------------------
library(plotly)
source('~/marketVirEvol/code/general/gen_functions.R')

# 1) set the fixed parameters --------------------------------------------------
virs = seq(1, 1000, 1)
c3 = 1e-03
b = 0
nat_mort = 1 / 365
gamma = 1 / 5
sigma = 1 / 5
N_m = 1000
p = 0
m_f = 0.1 / 120
N_f = 1e6
prev_f = 0.12
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.402
R_f = seroprev_f * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f
DFE_markets = (S_f * m_f) / (nat_mort + (1 / 5.5))

# Environmental res params
epsilon = 1
psi = 1 / 5
phi = 1

# 2) set the ranges for the parameters to vary ---------------------------------
c1_range = c(1 / 2300000, 1 / 230000, 1 / 23000, 1 / 2300, seq(1 / 230, 1 / 23, by=1/200), 1 / 22)
c2_range = seq(0.1, 1, 0.15)
market_psi_clean = 1
market_m_m = 1 / 5.5
psi_clean_range = seq(1, 10, 1) # This is equivalent to kappa in the text
m_m_range = c(seq(1 / 365, 1 / 5.5, 0.025), 1/5.5)
length(c1_range) * length(c2_range) * length(psi_clean_range) * length(m_m_range)
if (!all(m_m_range == cummax(m_m_range)) | !all(c1_range == cummax(c1_range)) | 
    !all(c2_range == cummax(c2_range)) | !all(psi_clean_range == cummax(psi_clean_range))) {
  stop('Must be in increasing order for tests.')
}

# 3) Find differences between farm conditions and market conditions for each viable tradeoff curve
# Vector used to find all differences between farm conditions and market conditions
diff_virs <- list()
diff_virs_dex <- 1
# Loop used to get answer to m_m
exclude_beta_cnt <- 0
for (c2 in c2_range) {
  print(paste0(round((which(round(c2_range, 2) == round(c2, 2)) - 1) / length(c2_range) * 100, 2),'% done.'))
  for (c1 in c1_range) {
    test_beta <- c1 * max(virs)^c2 * DFE_markets
    R0s <- c()
    for (vir in virs) {
      mort <- (vir) * c3
      beta <- c1 * (vir)^c2
      lambda <- beta * phi
      R0_markets <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                                sigma=sigma, nat_mort=nat_mort, m=market_m_m, gamma=gamma, mort=mort, kappa=market_psi_clean, psi=psi, both=T)
      R0s <- c(R0s, R0_markets)
    }
    opt_vir <- virs[which(R0s == max(R0s))]
    full_condition <- test_beta <= 200 & test_beta >= 1 & max(R0s) <= 100 & max(R0s) >= 1
    if (full_condition) { # If full condition met, then do all comparisons against the market values
      final_df <- as.data.frame(matrix(nrow=length(m_m_range), ncol=length(psi_clean_range)))
      for (m_dex in 1:length(m_m_range)) {
        for (k_dex in 1:length(psi_clean_range)) {
          R0s_compare <- c()
          for (vir in virs) {
            mort <- (vir) * c3
            beta <- c1 * (vir)^c2
            lambda <- beta * phi
            R0_compare <- get_R0_shed_diff_migrate(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                                                   sigma=sigma, nat_mort=nat_mort, m=m_m_range[m_dex], m_I=m_m_range[m_dex]*0.5, gamma=gamma, mort=mort, kappa=psi_clean_range[k_dex], psi=psi, both=T)
            R0s_compare <- c(R0s_compare, R0_compare)
          }
          opt_vir_compare <- virs[which(R0s_compare == max(R0s_compare))]
          final_df[m_dex, k_dex] <- opt_vir - opt_vir_compare
        }
      }
      diff_virs[[diff_virs_dex]] <- final_df
      diff_virs_dex <- diff_virs_dex + 1
    } else {
      exclude_beta_cnt <- exclude_beta_cnt + 1
    }
  }
}
# Save objects
save(diff_virs, file = paste0("~/marketVirEvol/code_output/obj/dens_shed_diff_migrate_all_", phi, ".RData"))
load(paste0("~/marketVirEvol/code_output/obj/dens_shed_diff_migrate_all_", phi, ".RData"))
