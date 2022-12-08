surv1 <- read.csv('~/Desktop/Survey TMM2.csv')
surv2 <- read.csv('~/Desktop/SurveyTNR3.csv')
mean(c(surv1$Days.until.sold, surv2$Days.until.sold))
# So m_m should be 3.54 days. However, this is once it arrives at the seller. The time from middleman
# to seller may also adde to this?

mean(c(surv1$chickenforsaleweekly, surv2$chickenforsaleweekly))
# So if there are 1000 poultry, this is roughly equivalent to about 16 sellers
# in a single market, if we assume that as soon as one poultry is sold, another
# is replaced in that market. What might be interesting is if these sellers are open year round,
# vs. if they are only open certain days of the week. These population bottleneck events of the host
# population of poultry may be important for the model, possibly. Therefore, we are seeking to 
# find whether a single market, or ten interconnected markets may possibly provide enough
# circulation to select for higher virulence.

# Writing is more important than the math, because conceptually, things have to make sense as well.





# Question 2:
# Solved for the R0 of the model, and assume that the ESS will be the R0. But we need to 
# assume some relationship between virulence and transmission and mortality.
# Set the factor to convert virulence to mortality.
# In this case, this means that mort will be equal to 1 in the worst case scenarios
# of any transmission curve, and will not cause any mortality in the best case scenarios
c1 <- 0.01
c2 <- 1
c3 <- 0.001
virs <- seq(0, 2000, 1)
betas <- c1 * (virs)^c2
plot(virs, betas, type='l')
plot(virs * c3, betas, type='l')

# Global parameters
N_m <- 1000
p <- 0.2
sigma <- 1 / 5
gamma <- 1 / 10
m_m <- 1 / 3.5
m_fI = m_f = 0.1 / 120

# Get E_f and I_f
N_f <- 1e6
prev_f <- 0.2
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.3
R_f = 0.3 * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f

# Model for simulation later
# No migration of infectious or incubation poultry from farms for this simulation since we are
# interested in whether a single infection can sustain transmission in the markets. Thus, this model
# is solely used for verification of whether transmission can occur in the simplest market model.
mod_eqn <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS_m = -(beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m +m_f*S_f -m_m*S_m
    dE_m = (beta/((S_m + E_m + I_m + R_m)^p))*S_m*I_m -m_m*E_m -sigma*E_m
    dI_m = sigma*E_m -m_m*I_m -gamma*I_m -mort*I_m
    dR_m = gamma*I_m +m_f*R_f -m_m*R_m
    return(list(c(dS_m, dE_m, dI_m, dR_m)))})}

# Get R0s
# Observation: R0 does not depend on the migration of infectious or incubation poultry into the system
R0s <- c()
for (vir in virs) {
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  R0 <- ((beta / (((m_f * (S_f)) / m_m))^(p)) * ((m_f * S_f) / m_m) * sigma) / ((-m_m -sigma) * (-m_m -gamma -mort))
  R0s <- c(R0s, R0)
}

# Plot
plot(virs, R0s, type='l')
plot(virs * c3, R0s, type='l')

# Try to show for two cases of N_f (1% and 0.1% of the farm population), and for all relevant c1, c2, c3 values, and m_m values, optimal R0 virulence increases with turnover rate
# The N_f values will create different equilibriums of the sizes of the market patch because of the migration rates
new_rows <- list()
new_rows_dex <- 1
new_rows_slow <- list()
new_rows_slow_dex <- 1
new_rows_med_slow <- list()
new_rows_med_slow_dex <- 1
final_res <- c()
m_m_res <- c()
flat_res <- c()
too_high_beta_cnt <- 0
for (N_f in c(1e5, 1e6)) {
  prev_f <- 0.2
  I_f = E_f = N_f * (prev_f / 2)
  seroprev_f <- 0.3
  R_f = 0.3 * (E_f + (1 - prev_f) * N_f)
  S_f = ((1 - prev_f) * N_f) - R_f
  for (c1 in seq(0.01, 2, 0.2)) {
    print(c1)
    for (c2 in seq(0.1, 1, 0.1)) {
      test_beta <- c1 * max(virs) ^c2
      if (test_beta <= 100 & test_beta >= 1) {
        for (c3 in c(0.001, 0.1, 1)) {
          # This innermost loop loops through the turnover rates from slow to fast
          # opt_virs will collect the optimal virulence strategies from slow to fast
          opt_virs <- c()
          res_R0s <- matrix(nrow=length(seq(1 / 365, 1 / 3.5, 0.1)), ncol=length(virs))
          res_R0s_dex <- 1
          for (m_m in seq(1 / 365, 1 / 3.5, 0.1)) {
            R0s <- c()
            # Go through all virulence strategies to find the strategy that is optimal
            for (vir in virs) {
              mort <- (vir) * c3
              beta <- c1 * (vir)^c2
              R0 <- ((beta / (((m_f * (S_f)) / m_m))^(p)) * ((m_f * S_f) / m_m) * sigma) / ((-m_m -sigma) * (-m_m -gamma -mort))
              R0s <- c(R0s, R0)
            }
            res_R0s[res_R0s_dex,] <- R0s
            res_R0s_dex <- res_R0s_dex + 1
            opt_vir <- virs[which(R0s == max(R0s))]
            opt_virs <- c(opt_virs, opt_vir)
            
            # Check that it will always monotonically decrease both before and after the optimal virulence
            # strategy (which also checks that there is only a single optimum)
            aft_opt <- R0s[which(R0s == max(R0s)):length(R0s)]
            if (!all(aft_opt == cummin(aft_opt))) {
              stop('Does not monotonically decrease after optimum')
            }
            bef_opt <- R0s[1:(which(R0s == max(R0s)))]
            if (!all(bef_opt == cummax(bef_opt))) {
              stop('Does not monotonically increase before optimum')
            }
            
            # If this is the fast turnover parameter, then run a simulation of the optimal
            # virulence strategy to see if the optimal virulence strategy will go extinct in
            # the market patch
            if (m_m == max(seq(1 / 365, 1 / 3.5, 0.1))) {
              # Run simulation with fast turnover
              mort <- opt_vir * c3
              beta <- c1 * (opt_vir)^c2
              parameters <- c(m_m = 1 / (3.5), m_fI=m_fI, gamma = gamma, m_f = m_f, beta = beta,
                              S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f,
                              p = p, mort = mort)
              init <- c(S_m = 999,
                        E_m = 0,
                        I_m = 1,
                        R_m = 0)
              times <- seq(1, 730, 1)
              out <- ode(y=init, times=times, mod_eqn, parms=parameters)
              out.df <- as.data.frame(out)
              extinct_fast <- F
              if (out.df$I_m[nrow(out.df)] < 1) {
                extinct_fast <- T
              }
              new_row <- data.frame(matrix(c(c1, c2, c3, opt_vir, max(R0s), extinct_fast), nrow=1, ncol=6))
              new_rows[[new_rows_dex]] <- new_row
              new_rows_dex <- new_rows_dex + 1
            } else if (m_m == min(seq(1 / 365, 1 / 3.5, 0.1))) { # If this is instead the slow turnover rate, save the maximum R0 achieved by the optimal virulence strategy
              # Run simulation with slow turnover
              mort <- opt_vir * c3
              beta <- c1 * (opt_vir)^c2
              parameters <- c(m_m = (1 / 365), m_fI=m_fI, gamma = gamma, m_f = m_f, beta = beta,
                              S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, sigma=sigma,
                              p = p, mort = mort)
              init <- c(S_m = 999,
                        E_m = 0,
                        I_m = 1,
                        R_m = 0)
              times <- seq(1, 730, 1)
              out <- ode(y=init, times=times, mod_eqn, parms=parameters)
              out.df <- as.data.frame(out)
              extinct_slow <- F
              if (out.df$I_m[nrow(out.df)] < 1) {
                extinct_slow <- T
              }
              new_row_slow <- data.frame(matrix(c(c1, c2, c3, opt_vir, max(R0s), extinct_slow), nrow=1, ncol=6))
              new_rows_slow[[new_rows_slow_dex]] <- new_row_slow
              new_rows_slow_dex <- new_rows_slow_dex + 1
            } else if (m_m == seq(1 / 365, 1 / 3.5, 0.1)[2]) {
              # Run simulation with medium slow turnover
              mort <- opt_vir * c3
              beta <- c1 * (opt_vir)^c2
              parameters <- c(m_m = (seq(1 / 365, 1 / 3.5, 0.1)[2]), m_fI=m_fI, gamma = gamma, m_f = m_f, beta = beta,
                              S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, sigma=sigma,
                              p = p, mort = mort)
              init <- c(S_m = 999,
                        E_m = 0,
                        I_m = 1,
                        R_m = 0)
              times <- seq(1, 730, 1)
              out <- ode(y=init, times=times, mod_eqn, parms=parameters)
              out.df <- as.data.frame(out)
              extinct_med_slow <- F
              if (out.df$I_m[nrow(out.df)] < 1) {
                extinct_med_slow <- T
              }
              new_row_med_slow <- data.frame(matrix(c(c1, c2, c3, opt_vir, max(R0s), extinct_med_slow), nrow=1, ncol=6))
              new_rows_med_slow[[new_rows_med_slow_dex]] <- new_row_med_slow
              new_rows_med_slow_dex <- new_rows_med_slow_dex + 1
            }
          }
          # Check that for all virulence strategies, R0 decreases as turnover rate increases
          m_m_res_final <- T
          for (col_dex in 1:ncol(res_R0s)) {
            if (!all(res_R0s[,col_dex] == cummin(res_R0s[,col_dex]))) {
              m_m_res_final <- F
            } 
          }
          m_m_res <- c(m_m_res, m_m_res_final)
          
          # Check that for all virulence strategies, the difference between the virulence
          # at max vir (if this is not the optimal) and the optimal virulence is smaller
          # when turnover is faster than when turnover is slower
          to_check <- c()
          for (m_m_dex in 1:length(opt_virs)) {
            if (opt_virs[m_m_dex] != max(virs)) {
              to_subtract_from <- res_R0s[m_m_dex, opt_virs[m_m_dex] + 1]
              to_subtract <- res_R0s[m_m_dex, ncol(res_R0s)]
              to_check <- c(to_check, to_subtract_from - to_subtract)
            }
          }
          if (length(to_check) >= 2) {
            flat_res <- c(flat_res, all(to_check == cummin(to_check)))
          }
          
          # Check if it's true that as turnover rate gets faster, optimal virulence strategy increases
          # for this set of c1, c2, and c3
          final_res <- c(final_res, all(opt_virs == cummax(opt_virs)))
        }
      } else {
        too_high_beta_cnt <- too_high_beta_cnt + 1
      }
    }
  }
}
extincts.df_fast <- do.call(rbind, new_rows)
colnames(extincts.df_fast) <- c('c1', 'c2', 'c3', 'opt_vir', 'R0', 'extinct_fast')
extincts.df_slow <- do.call(rbind, new_rows_slow)
colnames(extincts.df_slow) <- c('c1', 'c2', 'c3', 'opt_vir', 'R0', 'extinct_slow')
extincts.df_med_slow <- do.call(rbind, new_rows_med_slow)
colnames(extincts.df_med_slow) <- c('c1', 'c2', 'c3', 'opt_vir', 'R0', 'extinct_med_slow')

# Percent of discarded parameter sets
too_high_beta_cnt / (length(seq(0.01, 2, 0.2)) * length(seq(0.1, 1, 0.1)))

# If this is true, then the intuition is correct and the optimal R0 also increases 
# with increasing turnover rate. 
# I've justed things and it is true for the ranges I've tested.
all(final_res)

# Check if first part is correct:
# This is true, so it does seem like R0 falls as m_m increases
all(m_m_res)

# If the following is true, then it does get flatter as m_m increases
all(flat_res)

# If wanting to plot single trajectories but just changing c1, c2, and c3
c1 <- 1.0
c2 <- 0.4
c3 <- 1e-04
opt_vir <- 50
m_m <- 1 / 3.5
mort <- opt_vir * c3
beta <- c1 * (opt_vir)^c2
parameters <- c(m_m = m_m, m_fI=m_fI, gamma = gamma, m_f = m_f, beta = beta,
                S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f,
                p = p, mort = mort)
init <- c(S_m = 999,
          E_m = 0,
          I_m = 0,
          R_m = 0)
times <- seq(1, 730, 1)
out <- ode(y=init, times=times, mod_eqn, parms=parameters)
out.df <- as.data.frame(out)
plot(out.df$time, out.df$S_m, type='l', col='blue', ylim=c(0, max(out.df$S_m)))
lines(out.df$time, out.df$E_m, type='l', col='green')
lines(out.df$time, out.df$I_m, type='l', col='red')
lines(out.df$time, out.df$R_m, type='l', col='black')

# Histogram of the max_R0s across all parameters when turnover is fast vs. slow
max_R0s_fasts <- extincts.df_fast$R0[which(extincts.df_fast$c3 == 0.001)]
max_R0s_slows <- extincts.df_slow$R0[which(extincts.df_slow$c3 == 0.001)]
b <- min(c(max_R0s_fasts, max_R0s_slows)) - 0.00001
e <- max(c(max_R0s_fasts, max_R0s_slows))
ax <- pretty(b:e, n = 25)
hist_fast <- hist(max_R0s_fasts, breaks = ax, plot = FALSE) 
hist_slow <- hist(max_R0s_slows, breaks = ax, plot = FALSE) 
lt.blue <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
lt.pink <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
plot(hist_fast, col = lt.pink, main='', xlab='Maximum R0s')
plot(hist_slow, col = lt.blue, add=T)

# Histogram of increases
hist(extincts.df_fast$opt_vir[which(extincts.df_fast$c3 == 0.001)] - extincts.df_slow$opt_vir[which(extincts.df_slow$c3 == 0.001)], breaks=20)

# So we've shown two things: (1) the R0 will go lower for all virulence strategies when
# turnover rate increases (2) for the parameter spaces we study of concave transmission curves,
# the optimal virulence strategy always rises with increasing R0.

# In many instances, the consequence of (1) is that transmission can not really be sustained
# solely in the market patch because the R0 is less than 1, or very close to 1. This is perhaps 
# saying that while the contact rate is very high in markets, and transmission might be greater, 
# the high turnover rate also prevents infectious poultry from coming into contact with other
# poultry for transmission to happen because it happens so fast.

# Part 3: next question: can markets act as a relevant source of selection for high virulence
# poultry when we look at actual transmission in the models.

# Hypothesis: No. We show through our model simulating many relevant tradeoff curves that the
# optimal virulence strategy cannot sustain transmission. If transmission is slowed however,
# then it could have, and adding births could also disrupt it.

# The following shows that more of the optimal strategies go extinct with fast turnover
# than with slow turnover. Thus this means that most virulence strategies generally, will
# not be able to subsist completely on the market patch, it needs to find another population
# and cannot act as a single source. We need to verify and compare this to a model with
# births I think, to act as the farm patch. But still, around 25% of the parameters tested,
# it is possible for optimal virulence strategies to circulate solely in the market patch.
# Thus, in some cases, maybe it can really act as a source. Should gain some intuition for
# what these tradeoff curves look like and if they're even realistic. Maybe should do a screening
# process for those tradeoff curves where the beta is greater than 100?

# Actually I just increased the range, and if c1 is allowed to go up to 10, then actually the majority
# do not go extinct. Can be the case that it needs a very dramatic push up to see the effects. So if
# beta is high even for low virulent strategies, then this is no good, should discard these as well.
# Wouldn't make sense to say that even with low virulence it can infect 10 poultry a day.
length(which(extincts.df_fast$extinct_fast == 1)) / nrow(extincts.df_fast)
length(which(extincts.df_slow$extinct_slow == 1)) / nrow(extincts.df_slow)
length(which(extincts.df_med_slow$extinct_med_slow == 1)) / nrow(extincts.df_med_slow)

# Just want to check that as m_m increases, the following expression increases
m_ms = seq(1/365, 1/3.5, by=0.001)
res_m_ms = c()
for (m_m in m_ms) {
  res_m_ms = c(res_m_ms, m_m^2 + gamma*m_m + mort * m_m + sigma * m_m + sigma * gamma + sigma * mort)
}
plot(m_ms, res_m_ms, type='l')

