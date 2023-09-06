library(deSolve)
library(foreach)
library(doParallel)
source('~/marketVirEvol/code/general/gen_functions.R')

mod_eqn_nob <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS = -(beta_res/((S + E + I + R)^p))*S*I -epsilon*(beta_res/((S + E + I + R)^p))*S*H -nat_mort*S -m*S +m_f*S_f -(beta_inv/((S + E2 + I2 + R2)^p))*S*I2 -epsilon*(beta_inv/((S + E2 + I2 + R2)^p))*S*H2
    
    dE = (beta_res/((S + E + I + R)^p))*S*I +epsilon*(beta_res/((S + E + I + R)^p))*S*H -sigma*E -nat_mort*E -m*E
    dI = sigma*E -m*I -gamma*I -mort_res*I -nat_mort*I
    dR = gamma*I -m*R -nat_mort*R
    dH = lambda_res*I -psi*kappa*H
    
    dE2 = (beta_inv/((S + E2 + I2 + R2)^p))*S*I2 +epsilon*(beta_inv/((S + E2 + I2 + R2)^p))*S*H2 -sigma*E2 -nat_mort*E2 -m*E2
    dI2 = sigma*E2 -m*I2 -gamma*I2 -mort_inv*I2 -nat_mort*I2
    dR2 = gamma*I2 -m*R2 -nat_mort*R2
    dH2 = lambda_inv*I2 -psi*kappa*H2
    return(list(c(dS, dE, dI, dR, dH, dE2, dI2, dR2, dH2)))})}

c1 = 1 / 2300
c2 = 0.4
c3 = 1 / 1000
b = 0
p = 0
m_f = 0.1 / 120
N_f = 1e6
prev_f = 0.12
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.402
R_f = seroprev_f * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f
epsilon = 1
phi = 1
psi = 1 / 5
m = 1 / 5.5
gamma = 1 / 5
mu = 1 / 365
sigma = 1 / 5

# Lower cleaning
alphas = seq(1, 1000, 1)
virs = alphas
R0s <- c()
kappa = 1
for (vir in virs) {
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  lambda <- beta * phi
  R0 <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                    sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort, kappa=kappa, psi=psi, both=T)
  R0s <- c(R0s, R0)
}
opt_vir_low <- virs[which(R0s == max(R0s))]
opt_vir_low

# Higher cleaning
alphas = seq(1, 1000, 1)
virs = alphas
R0s <- c()
kappa = 9
for (vir in virs) {
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  lambda <- beta * phi
  R0 <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                    sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort, kappa=kappa, psi=psi, both=T)
  R0s <- c(R0s, R0)
}
opt_vir_high <- virs[which(R0s == max(R0s))]
opt_vir_high

if (opt_vir_low <= opt_vir_high) {
  stop('Optimal virulence should be lower when cleaning is higher, and the 
        virulence parameters should be chosen such that it is strictly lower.')
}

# The opt_vir should invade all other strategies
test_invade <- function(alpha, opt_vir_input, kappa) {
  if (alpha != opt_vir_input) {
    # Run resident
    mort_res <- (alpha) * c3
    beta_res <- c1 * (alpha)^c2
    lambda_res <- beta_res * phi
    R0_res <- get_R0_shed(beta=beta_res, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda_res,
                sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort_res, kappa=kappa, psi=psi, both=T)
    R0_res
    parameters <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                    S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                    p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi, kappa = kappa,
                    beta_inv = 0, lambda_inv = 0, mort_inv = 0)
    init <- c(S = 9999, 
              E = 0, 
              I = 1, 
              R = 0, 
              H = 0,
              E2 = 0, 
              I2 = 0, 
              R2 = 0, 
              H2 = 0)
    times <- seq(1, 1000, 1)
    out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
    out.df <- as.data.frame(out)
    
    # Do quick testing
    if (out.df[nrow(out.df),7] > 0 | out.df[nrow(out.df),8] > 0 | out.df[nrow(out.df),9] > 0 | out.df[nrow(out.df),10] > 0) {
      stop('Error in simulation.')
    }
    
    # Run invasion
    mort_inv <- (opt_vir_input) * c3
    beta_inv <- c1 * (opt_vir_input)^c2
    lambda_inv <- beta_inv * phi
    R0_inv <- get_R0_shed(beta=beta_inv, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda_inv,
                sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort_inv, kappa=kappa, psi=psi, both=T)
    R0_inv
    parameters2 <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                    S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                    p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi, kappa = kappa,
                    beta_inv = beta_inv, lambda_inv = lambda_inv, mort_inv = mort_inv)
    init2 <- c(S = out.df[nrow(out.df),2],
               E = out.df[nrow(out.df),3],
               I = out.df[nrow(out.df),4],
               R = out.df[nrow(out.df),5],
               H = out.df[nrow(out.df),6],
               E2 = 0,
               I2 = 1,
               R2 = 0,
               H2 = 0)
    times2 <- seq(1, 1e8, 100000)
    out2 <- ode(y=init2, times=times2, mod_eqn_nob, parms=parameters2)
    out.df2 <- as.data.frame(out2)
    #plot(times2, out.df2$I, type='l', ylim=c(0, max(out.df2$I, out.df2$I2)))
    #lines(times2, out.df2$I2, type='l', col='red')
    
    return(((out.df2[nrow(out.df2),3] + out.df2[nrow(out.df2),4]) < 1) & ((out.df2[nrow(out.df2),7] + out.df2[nrow(out.df2),8]) > 2))
  } else {
    return(TRUE)
  }
}

# The opt_vir should resist invasion from all other strategies
test_ES <- function(alpha, opt_vir_input, kappa) {
  if (alpha != opt_vir_input) {
    # Run resident
    mort_res <- (opt_vir_input) * c3
    beta_res <- c1 * (opt_vir_input)^c2
    lambda_res <- beta_res * phi
    get_R0_shed(beta=beta_res, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda_res,
                sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort_res, kappa=kappa, psi=psi, both=T)
    parameters <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                    S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                    p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi, kappa = kappa,
                    beta_inv = 0, lambda_inv = 0, mort_inv = 0)
    init <- c(S = 999, E = 0, I = 1, R = 0, H = 0,
              E2 = 0, I2 = 0, R2 = 0, H2 = 0)
    times <- seq(1, 1000, 1)
    out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
    out.df <- as.data.frame(out)
    
    # Do quick testing
    if (out.df[nrow(out.df),7] > 0 | out.df[nrow(out.df),8] > 0 | out.df[nrow(out.df),9] > 0 | out.df[nrow(out.df),10] > 0) {
      stop('Error in simulation.')
    }
    
    # Run invasion
    mort_inv <- (alpha) * c3
    beta_inv <- c1 * (alpha)^c2
    lambda_inv <- beta_inv * phi
    get_R0_shed(beta=beta_inv, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda_inv,
                sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort_inv, kappa=kappa, psi=psi, both=T)
    parameters2 <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                    S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                    p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi, kappa = kappa,
                    beta_inv = beta_inv, lambda_inv = lambda_inv, mort_inv = mort_inv)
    init2 <- c(S = out.df[nrow(out.df),2],
               E = out.df[nrow(out.df),3],
               I = out.df[nrow(out.df),4],
               R = out.df[nrow(out.df),5],
               H = out.df[nrow(out.df),6],
               E2 = 0,
               I2 = 1,
               R2 = 0,
               H2 = 0)
    times2 <- seq(1, 1e8, 100000)
    out2 <- ode(y=init2, times=times2, mod_eqn_nob, parms=parameters2)
    out.df2 <- as.data.frame(out2)
    
    return(((out.df2[nrow(out.df2),3] + out.df2[nrow(out.df2),4]) > 2) & ((out.df2[nrow(out.df2),7] + out.df2[nrow(out.df2),8]) < 1))
  } else {
    return(TRUE)
  }
}

# opt_vir invading -------------------------------------------------------------
# Run in parallel, invasion (low cleaning)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=seq(1, 1000, by=1), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(alpha=i, opt_vir_input=opt_vir_low, kappa = 1)
  tempMatrix
}
stopCluster(cl)
res1 <- seq(1, 1000, by=1)[which(!finalMatrix)]

# Run in parallel, invasion (high cleaning)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix2 <- foreach(i=seq(1, 1000, by=1), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(alpha=i, opt_vir_input=opt_vir_high, kappa = 9)
  tempMatrix
}
stopCluster(cl)
res2 <- seq(1, 1000, by=1)[which(!finalMatrix2)]

# Attempted invasion of opt_vir ------------------------------------------------
# Run in parallel, invasion (low cleaning)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix3 <- foreach(i=seq(1, 1000, by=1), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_ES(alpha=i, opt_vir_input=opt_vir_low, kappa = 1)
  tempMatrix
}
stopCluster(cl)
res3 <- seq(1, 1000, by=1)[which(!finalMatrix3)]

# Run in parallel, invasion (high cleaning)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix4 <- foreach(i=seq(1, 1000, by=1), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_ES(alpha=i, opt_vir_input=opt_vir_high, kappa = 9)
  tempMatrix
}
stopCluster(cl)
res4 <- seq(1, 1000, by=1)[which(!finalMatrix4)]

if (length(res1) != 0 & length(res2) != 0 & length(res3) != 0 & length(res4) != 0) {
  stop('Error.')
} else {
  print('Invasion analysis test successful.')
}
