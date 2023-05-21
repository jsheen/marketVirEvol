# ------------------------------------------------------------------------------
# @description: this is used to evaluate the second part of the virulence 
#               evolution of markets paper. It is divided into two parts: 
#               - (a) all questions besides extinction question, which requires simulation
#               - (b) extinction question, which requires simulation
# ------------------------------------------------------------------------------
# 0) libraries and sources -----------------------------------------------------
library(deSolve)
library(plotly)
source('~/marketVirEvol/code/general/gen_functions.R')

# 1) set models to be used for the extinction analysis -------------------------
# -------------------------------------------------------------------------------
# This model excludes the birth term, and will be used instead of the previous model due to a division by 0 error when b=0
mod_eqn_nob <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS_m = -(beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m -epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -nat_mort*S_m -m_m*S_m +m_f*S_f
    dE_m = (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m +epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m -sigma*E_m -nat_mort*E_m -m_m*E_m
    dI_m = sigma*E_m -m_m*I_m -gamma*I_m -mort*I_m -nat_mort*I_m
    dR_m = gamma*I_m -m_m*R_m -nat_mort*R_m
    dH_m = sigma*E_m -psi*kappa*H_m
    return(list(c(dS_m, dE_m, dI_m, dR_m, dH_m)))})}

# This is the full model that includes the birth term
mod_eqn_b <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS_m = b * (S_m + E_m + I_m + R_m) -  
      (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m - nat_mort * S_m + m_f * S_f - m_m * S_m -epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m
    dE_m = (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m -sigma*E_m -nat_mort*E_m -m_m*E_m +epsilon*(beta/((S_m + E_m + I_m + R_m)^p))*S_m*H_m
    dI_m = sigma*E_m -m_m*I_m -gamma*I_m -mort*I_m -nat_mort*I_m
    dR_m = gamma*I_m -m_m*R_m -nat_mort*R_m
    dH_m = sigma*E_m -psi*kappa*H_m
    return(list(c(dS_m, dE_m, dI_m, dR_m, dH_m)))})}

# 2) set the fixed parameters --------------------------------------------------
virs = seq(0, 1000, 1)
c3 = 1e-03
b = 0
nat_mort = 1 / 365
gamma = 1 / 5
sigma = 1 / 5
epsilon = 1
N_m = 1000
p = 0
m_f = 0.1 / 120
N_f = 1e6
prev_f = 0.12
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.402
R_f = seroprev_f * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f

# 3) set the ranges for the parameters to vary ---------------------------------
c1_range = c(1 / 2300000, 1 / 230000, 1 / 23000, 1 / 2300, seq(1 / 230, 1 / 23, by=1/200), 1 / 22)
c2_range = seq(0.1, 1, 0.15)
psi_clean_range = seq(1, 10, 2)
m_m_range = seq(1 / 365, 1 / 5.5, 0.025)
length(c1_range) * length(c2_range) * length(psi_clean_range) * length(m_m_range)
if (!all(m_m_range == cummax(m_m_range)) | !all(c1_range == cummax(c1_range)) | 
    !all(c2_range == cummax(c2_range)) | !all(psi_clean_range == cummax(psi_clean_range))) {
  stop('Must be in increasing order for the rest of the tests to make sense.')
}

# 5) run loop in order to see what percentage go extinct even if R0 >= 1 -------
extincts_full <- c()
extincts_full_less_one <- c()
final_Ns <- c()
equil_noDis_Ns <- c()
R0_less_one <- 0
should_be_extinct <- 0
for (c2 in c2_range) {
  print(paste0(round((which(round(c2_range, 2) == round(c2, 2)) - 1) / length(c2_range) * 100, 2),'% done.'))
  for (c1 in c1_range) {
    test_beta <- c1 * max(virs)^c2
    if (test_beta <= 100 & test_beta >= 1) {
      for (psi_clean in psi_clean_range) {
        for (m_m in m_m_range) {
          R0s <- c()
          for (vir in virs) {
            mort <- (vir) * c3
            beta <- c1 * (vir)^c2
            psi <-  (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
            R0 <- get_R0(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, sigma=sigma, nat_mort=nat_mort, m_m=m_m, gamma=gamma, mort=mort, psi=psi, kappa=psi_clean, both=F)
            R0s <- c(R0s, R0)
          }
          # What was the optimal virulence strategy for this m_m?
          opt_vir <- virs[which(R0s == max(R0s))]
          # If this optimal virulence strategy leads to R0 >= 1
          if (max(R0s) >= 1) {
            mort <- opt_vir * c3
            beta <- c1 * (opt_vir)^c2
            psi <-  (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
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
            extinct <- F
            if (out.df$I_m[nrow(out.df)] < 1) {
              extinct <- T
              # print(paste0('c1: ', c1))
              # print(paste0('c2: ', c2))
              # print(paste0('m_m: ', m_m))
              # print(paste0('psi_clean: ', psi_clean))
              # print(max(R0s))
              # print(opt_vir)
            }
            extincts_full <- c(extincts_full, paste0(c1, ',', c2, ',', m_m, ',', psi_clean, ',', extinct, ',', opt_vir))
            final_Ns <- c(final_Ns, sum(out.df[nrow(out.df),2:5]))
            equil_noDis_Ns <- c(equil_noDis_Ns, (m_f * N_f) / (nat_mort + m_m))
            if (((m_f * N_f) / (nat_mort + m_m)) - sum(out.df[nrow(out.df),2:5]) < 0) {
              stop('Error in equilibrium N.')
            }
          } else {
            R0_less_one <- R0_less_one + 1
            
            # Check if it goes extinct if it is less than one
            mort <- opt_vir * c3
            beta <- c1 * (opt_vir)^c2
            psi <- (1 / ((1 / (gamma + nat_mort + m_m + mort)) + (1 / gamma)))
            parameters <- c(m_m = m_m, gamma = gamma, m_f = m_f, beta = beta,
                            S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = nat_mort,
                            p = p, mort = mort, sigma = sigma, epsilon = epsilon, psi=psi, kappa=psi_clean)
            init <- c(S_m = 999,
                      E_m = 0,
                      I_m = 1,
                      R_m = 0,
                      H_m = 0)
            times <- seq(1, 730, 1)
            out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
            out.df <- as.data.frame(out)
            if (out.df$I_m[nrow(out.df)] < 1) {
              should_be_extinct <- should_be_extinct + 1
              extincts_full_less_one <- c(extincts_full_less_one, paste0(c1, ',', c2, ',', m_m, ',', psi_clean, ',', T, ',', opt_vir, ',', max(R0s)))
            } else {
              extincts_full_less_one <- c(extincts_full_less_one, paste0(c1, ',', c2, ',', m_m, ',', psi_clean, ',', F, ',', opt_vir, ',', max(R0s)))
            }
          }
        }
      }
    }
  }
}
# Save objects
save(extincts_full, extincts_full_less_one, final_Ns, equil_noDis_Ns, R0_less_one, should_be_extinct, file = "~/marketVirEvol/code_output/obj/extincts_dens_epsilon1.RData")
load("~/marketVirEvol/code_output/obj/extincts_dens_epsilon1.RData")

# Double checked that if R0 is less than 1, then it goes extinct
extincts_res <- sapply(extincts_full_less_one, function(x) as.logical(strsplit(x, ',')[[1]][5]))
length(which(extincts_res)) / length(extincts_res) * 100

# 6) answer questions from the loop --------------------------------------------
# What is the percentage of parameter combinations with R0 >= 1 that go extinct?:
extincts_res <- sapply(extincts_full, function(x) as.logical(strsplit(x, ',')[[1]][5]))
length(which(extincts_res)) / length(extincts_res) * 100

# Which parameter combinations go extinct of those with R0 >= 1?
if ((length(which(extincts_res)) / length(extincts_res) * 100) != 0) {
  # First create dataframe for the extinct parameters
  subset_extinct <- which(extincts_res)
  extincts_c1 <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][1]))[subset_extinct]
  extincts_c2 <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][2]))[subset_extinct]
  extincts_m_m <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][3]))[subset_extinct]
  extincts_psi_clean <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][4]))[subset_extinct]
  extincts_opt_vir <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][6]))[subset_extinct]
  orig_extincts_df <- data.frame(matrix(c(extincts_c1, extincts_c2, extincts_psi_clean, extincts_m_m, extincts_opt_vir), ncol=5, nrow=length(extincts_c1), byrow=F))
  colnames(orig_extincts_df) <- c('c1', 'c2', 'psi_clean', 'm_m', 'opt_vir')
  new_extincts_ls <- list()
  new_extincts_ls_dex <- 1
  for (input_c1 in c1_range) {
    for (input_c2 in c2_range) {
      for (input_m_m in m_m_range) {
        test_if_exists <- which(round(orig_extincts_df$c1, 4) == round(input_c1, 4) & round(orig_extincts_df$c2, 2) == round(input_c2, 2) & round(orig_extincts_df$m_m, 10) == round(input_m_m, 10))
        if (length(test_if_exists) > 0) {
          sub <- orig_extincts_df[test_if_exists,]
          ave_psi_clean <- mean(sub$psi_clean)
          new_row <- data.frame(matrix(c(input_c1, input_c2, input_m_m, ave_psi_clean), nrow=1, ncol=4))
          new_extincts_ls[[new_extincts_ls_dex]] <- new_row
          new_extincts_ls_dex <- new_extincts_ls_dex + 1
        }
      }
    }
  }
  new_extincts_df <- do.call(rbind, new_extincts_ls)
  colnames(new_extincts_df) <- c('c1', 'c2', 'm_m', 'ave_psi_clean')
}

# Then create a dataframe for the not extinct parameters of those with R0 >= 1
subset_not_extinct <- which(!extincts_res)
not_extincts_c1 <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][1]))[subset_not_extinct]
not_extincts_c2 <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][2]))[subset_not_extinct]
not_extincts_m_m <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][3]))[subset_not_extinct]
not_extincts_psi_clean <- sapply(extincts_full, function(x) as.numeric(strsplit(x, ',')[[1]][4]))[subset_not_extinct]
orig_not_extincts_df <- data.frame(matrix(c(not_extincts_c1, not_extincts_c2, not_extincts_psi_clean, not_extincts_m_m), ncol=4, nrow=length(not_extincts_c1), byrow=F))
colnames(orig_not_extincts_df) <- c('c1', 'c2', 'psi_clean', 'm_m')
# Reduce to three dimensions instead of four
new_not_extincts_ls <- list()
new_not_extincts_ls_dex <- 1
for (input_c1 in c1_range) {
  for (input_c2 in c2_range) {
    for (input_m_m in m_m_range) {
      test_if_exists <- which(round(orig_not_extincts_df$c1, 4) == round(input_c1, 4) & round(orig_not_extincts_df$c2, 2) == round(input_c2, 2) & round(orig_not_extincts_df$m_m, 10) == round(input_m_m, 10))
      if (length(test_if_exists) > 0) {
        sub <- orig_not_extincts_df[test_if_exists,]
        ave_psi_clean <- mean(sub$psi_clean)
        new_row <- data.frame(matrix(c(input_c1, input_c2, input_m_m, ave_psi_clean), nrow=1, ncol=4))
        new_not_extincts_ls[[new_not_extincts_ls_dex]] <- new_row
        new_not_extincts_ls_dex <- new_not_extincts_ls_dex + 1
      }
    }
  }
}
new_not_extincts_df <- do.call(rbind, new_not_extincts_ls)
colnames(new_not_extincts_df) <- c('c1', 'c2', 'm_m', 'ave_psi_clean')

# First visualize extincts when R0 is greater than or equal to 1
fig <- plot_ly(data=new_extincts_df, x =~c1, y = ~c2, z = ~m_m, color=~ave_psi_clean)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'c1', range=c(min(c1_range),max(c1_range))), yaxis = list(title = 'c2', range=c(min(c2_range),max(c2_range))), zaxis = list(title = 'm_m', range=c(min(m_m_range),max(m_m_range)))))
fig <- fig %>% layout(scene = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                   yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                   zaxis = list(title = list(text='<b>m<sub>m</sub></b>', font=list(size=30)))))
fig

# # Then visualize not extincts just to double check that loop worked correctly
# So among not_extincts, these occur when c1 is low, and there is a not very flat curve, and can occur when m_m is a bit low but with a lower psi_clean, or m_m is a bit higher with a higher_psi_clean
fig <- plot_ly(data=new_not_extincts_df, x = ~c1, y = ~c2, z = ~m_m, color=~ave_psi_clean)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'c1', range=c(min(c1_range),max(c1_range))), yaxis = list(title = 'c2', range=c(min(c2_range),max(c2_range))), zaxis = list(title = 'm_m', range=c(min(m_m_range),max(m_m_range)))))
fig <- fig %>% layout(scene = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                   yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                   zaxis = list(title = list(text='<b>m<sub>m</sub></b>', font=list(size=30)))))

fig

# What percentage of combinations have R0 < 1 out of the parameters searched?
# Around 1.1% have R0 strictly less than 1
R0_less_one / (length(extincts_full) + R0_less_one) * 100

# What is the difference between the DFE equilibrium and the equilbrium with optimal virulence?
hist(equil_noDis_Ns - final_Ns, breaks=100)
