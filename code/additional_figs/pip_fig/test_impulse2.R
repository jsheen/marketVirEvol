library(deSolve)
library(foreach)
library(doParallel)
source('~/marketVirEvol/code/general/gen_functions.R')
start_time <- Sys.time()
# Model and parameters ---------------------------------------------------------
mod_eqn_nob <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dS = -(beta_res/((S + E + I + R)^p))*S*I -epsilon*(beta_res/((S + E + I + R)^p))*S*H -nat_mort*S -m*S +m_f*S_f -(beta_inv/((S + E2 + I2 + R2)^p))*S*I2 -epsilon*(beta_inv/((S + E2 + I2 + R2)^p))*S*H2
    
    dE = (beta_res/((S + E + I + R)^p))*S*I +epsilon*(beta_res/((S + E + I + R)^p))*S*H -sigma*E -nat_mort*E -m*E
    dI = sigma*E -m*I -gamma*I -mort_res*I -nat_mort*I
    dR = gamma*I -m*R -nat_mort*R
    dH = lambda_res*I -psi*H
    
    dE2 = (beta_inv/((S + E2 + I2 + R2)^p))*S*I2 +epsilon*(beta_inv/((S + E2 + I2 + R2)^p))*S*H2 -sigma*E2 -nat_mort*E2 -m*E2
    dI2 = sigma*E2 -m*I2 -gamma*I2 -mort_inv*I2 -nat_mort*I2
    dR2 = gamma*I2 -m*R2 -nat_mort*R2
    dH2 = lambda_inv*I2 -psi*H2
    return(list(c(dS, dE, dI, dR, dH, dE2, dI2, dR2, dH2)))})}

c1 = 1 / 2300
c2 = 0.45
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
phi = 10
psi = 1 / 5
m = 1 / 5.5
gamma = 1 / 5
mu = 1 / 365
sigma = 1 / 5

# Params for invasion analysis -------------------------------------------------
cohort_dur = 90
t_step = 1
n_gen = 1000
n_gen2 = 1000

# First invasion analysis when cleaning (kappa) is low -------------------------
test_invade <- function(alpha, inv_alpha, kappa) {
  if (alpha != inv_alpha) {
    # Run resident for n_gen number of generations
    mort_res <- (alpha) * c3
    beta_res <- c1 * (alpha)^c2
    lambda_res <- beta_res * phi
    parameters <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                    S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                    p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi,
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
    times <- c(seq(1, cohort_dur, t_step), cohort_dur)
    for (gen in 1:n_gen) {
      out <- ode(y=init, times=times, mod_eqn_nob, parms=parameters)
      out.df <- as.data.frame(out)
      init <- c(S = 9999,#out.df[nrow(out.df),2], 
                E = 0,#out.df[nrow(out.df),3],
                I = 0,#out.df[nrow(out.df),4], 
                R = 0,#out.df[nrow(out.df),5], 
                H = out.df[nrow(out.df),6] * (1 - kappa),
                E2 = 0,#out.df[nrow(out.df),7],
                I2 = 0,#out.df[nrow(out.df),8], 
                R2 = 0,#out.df[nrow(out.df),9], 
                H2 = out.df[nrow(out.df),10] * (1 - kappa))
    }
    
    # Run invasion
    # Do quick testing
    if (out.df[nrow(out.df),7] > 0 | out.df[nrow(out.df),8] > 0 | out.df[nrow(out.df),9] > 0 | out.df[nrow(out.df),10] > 0) {
      stop('Error in simulation.')
    }
    mort_inv <- (inv_alpha) * c3
    beta_inv <- c1 * (inv_alpha)^c2
    lambda_inv <- beta_inv * phi
    parameters2 <- c(m = m, gamma = gamma, m_f = m_f, beta_res = beta_res,
                     S_f = S_f, E_f = E_f, I_f = I_f, R_f = R_f, b = b, nat_mort = mu, lambda_res=lambda_res,
                     p = p, mort_res = mort_res, sigma = sigma, epsilon = epsilon, psi=psi,
                     beta_inv = beta_inv, lambda_inv = lambda_inv, mort_inv = mort_inv)
    init2 <- c(S = out.df[nrow(out.df),2],
               E = out.df[nrow(out.df),3],
               I = out.df[nrow(out.df),4], 
               R = out.df[nrow(out.df),5], 
               H = out.df[nrow(out.df),6],
               E2 = out.df[nrow(out.df),7],
               I2 = 1, 
               R2 = out.df[nrow(out.df),9],
               H2 = out.df[nrow(out.df),10])
    times2 <- c(seq(1, cohort_dur, t_step), cohort_dur)
    for (gen in 1:n_gen2) {
      out2 <- ode(y=init2, times=times2, mod_eqn_nob, parms=parameters2)
      out2.df <- as.data.frame(out2)
      init2 <- c(S = 9999,#out2.df[nrow(out2.df),2], 
                 E = 0,#out2.df[nrow(out2.df),3],
                 I = 0,#out2.df[nrow(out2.df),4], 
                 R = 0,#out2.df[nrow(out2.df),5], 
                 H = out2.df[nrow(out2.df),6] * (1 - kappa),
                 E2 = 0,#out2.df[nrow(out2.df),7],
                 I2 = 0,#out2.df[nrow(out2.df),8], 
                 R2 = 0,#out2.df[nrow(out2.df),9], 
                 H2 = out2.df[nrow(out2.df),10] * (1 - kappa))
    }
    
    # Last gen
    # times3 <- c(seq(1, cohort_dur * 10, (t_step * 10)), cohort_dur * 10)
    # out3 <- ode(y=init2, times=times3, mod_eqn_nob, parms=parameters2)
    # out3.df <- as.data.frame(out3)
    # out2.df <- out3.df
    
    if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) < -1) |
        ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) < -1)) {
      stop('Error.')
    }
    
    #plot(times2, out2.df$I, type='l', ylim=c(0, max(out2.df$I, out2.df$I2)))
    #lines(times2, out2.df$I2, type='l', col='red')
    if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) < 1) & ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) > 2)) {
      return(1)
    } else if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) >= 1) & ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) < 1)) {
      return(0)
    } else if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) >= 1) & (((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) > 1) & ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) <= 2))) {
      return(4)
    } else if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) >= 1) & ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) > 2))  {
      return(2)
    } else if (((out2.df[nrow(out2.df),3] + out2.df[nrow(out2.df),4]) < 1) & ((out2.df[nrow(out2.df),7] + out2.df[nrow(out2.df),8]) < 1)) {
      return(3)
    }
  } else {
    return(1)
  }
}

alphas = seq(1, 1000, 10)
all_comb = expand.grid(alphas, alphas)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=1:nrow(all_comb), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(all_comb[i,1], all_comb[i,2], kappa=0)
  tempMatrix
}
stopCluster(cl)
kapLowMat <- apply(matrix(finalMatrix, ncol=length(alphas), nrow=length(alphas), byrow=T), 2, rev)
kapLowMat
write.csv(kapLowMat, '~/Desktop/kapLowMat2.csv')

# Then invasion analysis when kappa is high ------------------------------------
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix2 <- foreach(i=1:nrow(all_comb), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(all_comb[i,1], all_comb[i,2], kappa=0.95)
  tempMatrix
}
stopCluster(cl)
kapHighMat <- apply(matrix(finalMatrix2, ncol=length(alphas), nrow=length(alphas), byrow=T), 2, rev)
kapHighMat
write.csv(kapHighMat, '~/Desktop/kapHighMat2.csv')

# Just to doublecheck by eye that matrix calculations are correct
toView <- c()
for (i in 1:nrow(all_comb)) {
  toView <- c(toView, paste0(all_comb$Var1[i], ',', all_comb$Var2[i]))
}
byEye <- apply(matrix(toView, ncol=length(alphas), nrow=length(alphas), byrow=T), 2, rev)
byEye

# There are a few cells that are neither 1 nor 0 in invasion analysis when kappa
# is high. Find these, and run for a longer number of generations --------------
kapLowMat <- read.csv('~/Desktop/stepSize10/kapLowMat2.csv')
kapLowMat <- kapLowMat[,2:ncol(kapLowMat)]
kapHighMat <- read.csv('~/Desktop/stepSize10/kapHighMat2.csv')
kapHighMat <- kapHighMat[,2:ncol(kapHighMat)]

n_gen = 1000
n_gen2 = 10000
which(kapLowMat != 1 & kapLowMat != 0, arr.ind=T)
lowCorrect <- which(kapLowMat != 1 & kapLowMat != 0, arr.ind=T)
lowComb <- c()
for (i in 1:nrow(lowCorrect)) {
  rowi <- lowCorrect[i,1]
  coli <- lowCorrect[i,2]
  lowComb <- c(lowComb, byEye[rowi, coli])
}
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
lowCorrected <- foreach(i=1:length(lowComb), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(as.numeric(strsplit(lowComb[i], ',')[[1]][1]), as.numeric(strsplit(lowComb[i], ',')[[1]][2]), kappa=0)
  tempMatrix
}
stopCluster(cl)
write.csv(lowCorrected, '~/Desktop/lowCorrected.csv', row.names=F)

which(kapHighMat != 1 & kapHighMat != 0, arr.ind=T)
highCorrect <- which(kapHighMat != 1 & kapHighMat != 0, arr.ind=T)
highComb <- c()
for (i in 1:nrow(highCorrect)) {
  rowi <- highCorrect[i,1]
  coli <- highCorrect[i,2]
  highComb <- c(highComb, byEye[rowi, coli])
}
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
highCorrected <- foreach(i=1:length(highComb), .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(as.numeric(strsplit(highComb[i], ',')[[1]][1]), as.numeric(strsplit(highComb[i], ',')[[1]][2]), kappa=0.95)
  tempMatrix
}
stopCluster(cl)
write.csv(highCorrected, '~/Desktop/highCorrected.csv', row.names=F)

# Fill in the new results ------------------------------------------------------
kapLowMat_corrected <- kapLowMat
for (i in 1:length(lowCorrected)) {
  kapLowMat_corrected[lowCorrect[i,1],lowCorrect[i,2]] <- lowCorrected[i]
  
}

kapHighMat_corrected <- kapHighMat
for (i in 1:length(highCorrected)) {
  kapHighMat_corrected[highCorrect[i,1],highCorrect[i,2]] <- highCorrected[i]
}

# One of them still did not resolve (index 7028):
# which(highCorrected == 4) : 135
# highComb[135] : '701,721'
n_gen = 1000
n_gen2 = 200000
test_invade(701, 721, kappa=0.95) #2
kapHighMat_corrected[highCorrect[135,1],highCorrect[135,2]] <- 2

write.csv(kapLowMat_corrected, '~/Desktop/stepSize10/kapLowMat2_corrected.csv', row.names = F)
write.csv(kapHighMat_corrected, '~/Desktop/stepSize10/kapHighMat2_corrected.csv', row.names = F)


# Find the optimal R0 using our model ------------------------------------------
source('~/marketVirEvol/code/general/gen_functions.R')
virs = seq(1, 1000, 1)
kappa_low = 1
kappa_high = 9
R0s_low <- c()
R0s_high <- c()
for (vir in virs) {
  mort <- (vir) * c3
  beta <- c1 * (vir)^c2
  lambda <- beta * phi
  R0_low <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                        sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort, kappa=kappa_low, psi=psi, both=T)
  R0_high <- get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                         sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=mort, kappa=kappa_high, psi=psi, both=T)
  R0s_low <- c(R0s_low, R0_low)
  R0s_high <- c(R0s_high, R0_high)
}
if (which(R0s_low == max(R0s_low)) > which(R0s_high == max(R0s_high))) {
  print('ESS lower as cleaning higher in original model.')
  print(which(R0s_low == max(R0s_low)))
  print(which(R0s_high == max(R0s_high)))
} else {
  print('Error.')
}

end_time <- Sys.time()
end_time - start_time

#rm(list = ls())

