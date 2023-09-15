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
phi = 10

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

# 3) set 1: m_m questions ------------------------------------------------------
# Q0. This vector is used to store the answer to whether there is a single optimum
opt_mm_res <- c()
# Q1. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, does the R0 for all virulence strategies decrease?
R0_mm_res <- c()
# Q2. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, is the difference between the max virulence and the optimal virulence smaller (that is, higher virulence strategies are more like lower virulence strategies)?
flat_mm_res <- c()
# Q3. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, does the optimal virulence strategy increase?
inc_mm_res <- c()
# Vector used to find all differences between slowest m_m and fastest m_m
diff_virs_m_m <- c()
# Loop used to get answer to m_m
exclude_beta_cnt <- 0
for (c2 in c2_range) {
  print(paste0(round((which(round(c2_range, 2) == round(c2, 2)) - 1) / length(c2_range) * 100, 2),'% done.'))
  for (c1 in c1_range) {
    test_beta <- c1 * max(virs)^c2 * DFE_markets
    full_condition <- test_beta <= 200 & test_beta >= 1
    if (full_condition) {
      # The following variables, opt_virs, and res_R0s, are used to answer store the answers to each question
      opt_virs <- c()
      res_R0s <- matrix(nrow=length(m_m_range), ncol=length(virs))
      res_R0s_dex <- 1
      # This innermost loop loops through each of the m_m for comparison
      for (m_m in m_m_range) {
        R0s <- c()
        for (vir in virs) {
          mort <- (vir) * c3
          beta <- c1 * (vir)^c2
          lambda <- beta * phi
          psi <- (1 / 5) * (vir) * c3
          R0 <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                            sigma=sigma, nat_mort=nat_mort, m=m_m, gamma=gamma, mort=mort, kappa=market_psi_clean, psi=psi, both=T)
          R0s <- c(R0s, R0)
        }
        # Store all virulence strategies for this m_m
        res_R0s[res_R0s_dex,] <- R0s
        res_R0s_dex <- res_R0s_dex + 1
        # What was the optimal virulence strategy for this m_m?
        opt_vir <- virs[which(R0s == max(R0s))]
        if (length(opt_vir) > 1) {
          stop('There should only be one maximum.')
        }
        opt_virs <- c(opt_virs, opt_vir)
        
        if (m_m == max(m_m_range)) {
          # Highest mm, optimal R0:
          opt_R0_max_mm = max(R0s)
        }
        else if (m_m == min(m_m_range)) {
          opt_R0_min_mm = max(R0s)
        }
        
        # Q0. Check if there is a single optimum
        aft_opt <- R0s[which(R0s == max(R0s)):length(R0s)]
        bef_opt <- R0s[1:(which(R0s == max(R0s)))]
        if (!all(aft_opt == cummin(aft_opt)) | !all(bef_opt == cummax(bef_opt))) {
          opt_mm_res <- c(opt_mm_res, F)
        } else {
          opt_mm_res <- c(opt_mm_res, T)
        }
      }
      
      # Q1. Check that for all virulence strategies, R0 decreases as turnover rate increases
      R0_mm_res_input <- T
      for (col_dex in 1:ncol(res_R0s)) {
        if (!all(res_R0s[,col_dex] == cummin(res_R0s[,col_dex]))) {
          R0_mm_res_input <- F
        }
      }
      if (!R0_mm_res_input) {
        stop('We are so confident this should be true, we will stop the entire loop if this is not true. The following parameter did not work.')
        print(paste0(c1, '; ', c2, '; ', c3))
      }
      R0_mm_res <- c(R0_mm_res, R0_mm_res_input)
      
      # Q2. Check gets flatter after optimum when increasing m_m
      to_check <- c()
      for (m_m_dex in 1:length(opt_virs)) {
        if (opt_virs[m_m_dex] != max(virs)) {
          # Get the R0 of the opt virulence strategy
          to_subtract_from <- res_R0s[m_m_dex, opt_virs[m_m_dex] + 1]
          # Get the R0 of the max virulence strategy
          to_subtract <- res_R0s[m_m_dex, ncol(res_R0s)]
          to_check <- c(to_check, to_subtract_from - to_subtract)
        }
      }
      if (length(to_check) >= 2) {
        flat_mm_res <- c(flat_mm_res, all(to_check == cummin(to_check)))
      }
      
      # Q3. Check whether optimal virulence increases when m_m increases
      inc_mm_res <- c(inc_mm_res, all(opt_virs == cummax(opt_virs)))
      
      # Find difference in optimal virulence strategies when m_m is slow vs. fast
      diff_virs_m_m <- c(diff_virs_m_m, paste0(c2, ',', c1, ',', market_psi_clean, ',', opt_virs[length(opt_virs)] - opt_virs[1], ',', opt_R0_max_mm, ',', opt_R0_min_mm))
      
    } else {
      exclude_beta_cnt <- exclude_beta_cnt + 1
    }
  }
}
# Save objects
save(opt_mm_res, R0_mm_res, flat_mm_res, inc_mm_res, diff_virs_m_m, file = paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_ref1_", lambda, "_", phi, ".RData"))
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_ref1_", lambda, "_", phi, ".RData"))

# 4) set 2: psi_clean_range questions ------------------------------------------
# Q0. This vector is used to store the answer to whether there is a single optimum
opt_psi_res <- c()
# Q1. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, does the R0 for all virulence strategies decrease?
R0_psi_res <- c()
# Q2. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, is the difference between the max virulence and the optimal virulence smaller (that is, higher virulence strategies are more like lower virulence strategies)?
flat_psi_res <- c()
# Q3. This vector is used to store the answer to: for a combination of c1, c2, and psi, as m_m increases, does the optimal virulence strategy increase?
inc_psi_res <- c()
# Loop used to get answer to psi
exclude_beta_cnt <- 0
# Vector used to find all differences between slowest psi and fastest psi
diff_virs_psi <- c()
for (c2 in c2_range) {
  print(paste0(round((which(round(c2_range, 2) == round(c2, 2)) - 1) / length(c2_range) * 100, 2),'% done.'))
  for (c1 in c1_range) {
    test_beta <- c1 * max(virs)^c2 * DFE_markets
    full_condition <- test_beta <= 200 & test_beta >= 1
    if (full_condition) {
      # The following variables, opt_virs, and res_R0s, are used to answer store the answers to each question
      opt_virs <- c()
      res_R0s <- matrix(nrow=length(psi_clean_range), ncol=length(virs))
      res_R0s_dex <- 1
      # This innermost loop loops through each of the psi for comparison
      for (psi_clean in psi_clean_range) {
        R0s <- c()
        for (vir in virs) {
          mort <- (vir) * c3
          beta <- c1 * (vir)^c2
          lambda <- beta * phi
          psi <- (1 / 5) * (vir) * c3
          R0 <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                            sigma=sigma, nat_mort=nat_mort, m=market_m_m, gamma=gamma, mort=mort, kappa=psi_clean, psi=psi, both=T)
          R0s <- c(R0s, R0)
        }
        
        # Store all virulence strategies for this m_m
        res_R0s[res_R0s_dex,] <- R0s
        res_R0s_dex <- res_R0s_dex + 1
        # What was the optimal virulence strategy for this m_m?
        opt_vir <- virs[which(R0s == max(R0s))]
        if (length(opt_vir) > 1) {
          stop('There should only be one maximum.')
        }
        opt_virs <- c(opt_virs, opt_vir)
        
        if (psi_clean == max(psi_clean_range)) {
          # Highest mm, optimal R0:
          opt_R0_max_kappa = max(R0s)
        }
        else if (psi_clean == min(psi_clean_range)) {
          opt_R0_min_kappa = max(R0s)
        }
        
        # Q0. Check if there is a single optimum
        aft_opt <- R0s[which(R0s == max(R0s)):length(R0s)]
        bef_opt <- R0s[1:(which(R0s == max(R0s)))]
        if (!all(aft_opt == cummin(aft_opt)) | !all(bef_opt == cummax(bef_opt))) {
          opt_psi_res <- c(opt_psi_res, F)
        } else {
          opt_psi_res <- c(opt_psi_res, T)
        }
      }
      
      
      # Q1. Check that for all virulence strategies, R0 decreases as psi increases
      R0_psi_res_input <- T
      for (col_dex in 1:ncol(res_R0s)) {
        if (!all(res_R0s[,col_dex] == cummin(res_R0s[,col_dex]))) {
          R0_psi_res_input <- F
        }
      }
      if (!R0_psi_res_input) {
        stop('We are so confident this should be true, we will stop the entire loop if this is not true. The following parameter did not work.')
        print(paste0(c1, '; ', c2, '; ', c3))
      }
      R0_psi_res <- c(R0_psi_res, R0_psi_res_input)
      
      # Q2. Check diff between R0 of opt vir and R0 of max vir gets sharper as psi increases
      to_check <- c()
      for (psi_dex in 1:length(opt_virs)) {
        if (opt_virs[psi_dex] != max(virs)) {
          # Get the R0 of the opt virulence strategy
          to_subtract_from <- res_R0s[psi_dex, opt_virs[psi_dex] + 1]
          # Get the R0 of the max virulence strategy
          to_subtract <- res_R0s[psi_dex, ncol(res_R0s)]
          to_check <- c(to_check, to_subtract_from - to_subtract)
        }
      }
      if (length(to_check) >= 2) {
        flat_psi_res <- c(flat_psi_res, all(to_check == cummax(to_check)))
      }
      
      # Q3. Check that opt vir decreases as psi increases
      inc_psi_res <- c(inc_psi_res, all(opt_virs == cummin(opt_virs)))
      
      # Find difference in optimal virulence strategies when m_m is slow vs. fast
      diff_virs_psi <- c(diff_virs_psi, paste0(c2, ',', c1, ',', market_m_m, ',', opt_virs[1] - opt_virs[length(opt_virs)], ',', opt_R0_max_kappa, ',', opt_R0_min_kappa))
      
    } else {
      exclude_beta_cnt <- exclude_beta_cnt + 1
    }
  }
}
# Save objects
save(opt_psi_res, R0_psi_res, flat_psi_res, inc_psi_res, diff_virs_psi, file = paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_ref1_", lambda, "_", phi, ".RData"))
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_ref1_", lambda, "_", phi, ".RData"))

# All plotting -----------------------------------------------------------------
max_R0_consider = 100

diff_virs_c2 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][1]))
diff_virs_c1 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][2]))
diff_virs_psi_clean <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][3]))
diff_virs_col <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
diff_virs_R02 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][6]))
constrain_R0 <- which(diff_virs_R0 <= max_R0_consider & diff_virs_R0 >= 1)# & diff_virs_R02 <= max_R0_consider & diff_virs_R02 >= 1)
fig1 <- plot_ly(x = diff_virs_c1[constrain_R0], y = diff_virs_c2[constrain_R0], z = diff_virs_psi_clean[constrain_R0], color=diff_virs_col[constrain_R0])
fig1 <- fig1 %>% add_markers()
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                     yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                     zaxis = list(title = list(text='<b>κ</b>', font=list(size=30)))))
fig1

# Plot 3D plot of differences
diff_virs_c2 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][1]))
diff_virs_c1 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][2]))
diff_virs_mm <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][3]))
diff_virs_col <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
diff_virs_R02 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][6]))
constrain_R0 <- which(diff_virs_R0 <= max_R0_consider & diff_virs_R0 >= 1)# & diff_virs_R02 <= max_R0_consider & diff_virs_R02 >= 1)
fig2 <- plot_ly(x = diff_virs_c1[constrain_R0], y = diff_virs_c2[constrain_R0], z = diff_virs_mm[constrain_R0], color=diff_virs_col[constrain_R0])
fig2 <- fig2 %>% add_markers()
fig2 <- fig2 %>% layout(scene2 = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                      yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                      zaxis = list(title = list(text='<b>m</b>', font=list(size=30)))))
fig2

# 6) Figure of the optimal R0s -------------------------------------------------
diff_virs_c2 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][1]))
diff_virs_c1 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][2]))
diff_virs_psi_clean <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][3]))
diff_virs_col <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
diff_virs_R02 <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][6]))
constrain_R0 <- which(diff_virs_R0 <= max_R0_consider & diff_virs_R0 >= 1 & diff_virs_R02 <= max_R0_consider & diff_virs_R02 >= 1)
fig3 <- plot_ly(x = diff_virs_c1[constrain_R0], y = diff_virs_c2[constrain_R0], z = diff_virs_psi_clean[constrain_R0], color=diff_virs_R0[constrain_R0])
fig3 <- fig3 %>% add_markers()
fig3 <- fig3 %>% layout(scene = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                     yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                     zaxis = list(title = list(text='<b>κ</b>', font=list(size=30)))))
fig3

diff_virs_c2 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][1]))
diff_virs_c1 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][2]))
diff_virs_mm <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][3]))
diff_virs_col <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
diff_virs_R02 <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][6]))
constrain_R0 <- which(diff_virs_R0 <= max_R0_consider & diff_virs_R0 >= 1 & diff_virs_R02 <= max_R0_consider & diff_virs_R02 >= 1)
fig4 <- plot_ly(x = diff_virs_c1[constrain_R0], y = diff_virs_c2[constrain_R0], z = diff_virs_mm[constrain_R0], color=diff_virs_R0[constrain_R0])
fig4 <- fig4 %>% add_markers()
fig4 <- fig4 %>% layout(scene2 = list(xaxis = list(title = list(text='<b>c<sub>1</sub></b>', font=list(size=30))), 
                                      yaxis = list(title = list(text='<b>c<sub>2</sub></b>', font=list(size=30))), 
                                      zaxis = list(title = list(text='<b>m</b>', font=list(size=30)))))
fig4





