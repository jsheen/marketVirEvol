# Libraries and source ---------------------------------------------------------
library(deSolve)
library(foreach)
library(doParallel)
library(pracma)
library(plot.matrix)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
source('~/virEvol/code/plot/plot_functions.R')

# General parameters -----------------------------------------------------------
# Time of simulation to first and second equilibria (days)
t_max_eq1 = t_max_eq2 = 4e3
# Population size
pop_size = 1e6
# Transition rate of infectiousness per chicken per day
sig = 1 / 5 
# Transition rate of recovery per chicken per day
gamm = 1 / 10
# Natural birth rate per chicken per day
b = 1 / 120
# Natural mortality rate per chicken per day
nat_mort = 1 / 365
# Mortality rate due to disease per chicken per day
mort = 1 / 4
# Virulence steps virus is allowed to take
vir_steps = seq(2, 100, 5)

# Model equation ---------------------------------------------------------------
eqn_mod4_freq <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet2*fS*fI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet1*fS*fVI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet2*fS*fVI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -v*fS +v_hat*fV -nat_mort*fS +m_mf*mS*(1-p_s) -m_fm*fS
    dfE1 = (fbet1*fS*fI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow +(fbet1*fS*fVI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -sig*fE1 -nat_mort*fE1 +m_mf*mE1*(1-p_s) -m_fm*fE1
    dfE2 = (fbet2*fS*fI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow +(fbet2*fS*fVI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -sig*fE2 -nat_mort*fE2 +m_mf*mE2*(1-p_s) -m_fm*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1 +m_mf*mI1*(1-p_s) -m_fm*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2 +m_mf*mI2*(1-p_s) -m_fm*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*fVI1 +gamm*fVI2 -nat_mort*fR +m_mf*mR*(1-p_s) -m_fm*fR
    dfV = -(fbet1*fV*fI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet2*fV*fI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet1*fV*fVI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -(fbet2*fV*fVI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow +v*fS -v_hat*fV -nat_mort*fV +m_mf*mV*(1-p_s) -m_fm_vax*fV
    dfVE1 = (fbet1*fV*fI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow +(fbet1*fV*fVI1)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -sig*fVE1 -nat_mort*fVE1 +m_mf*mVE1*(1-p_s) -m_fm_vax*fVE1
    dfVE2 = (fbet2*fV*fI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow +(fbet2*fV*fVI2)/(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)^N_pow -sig*fVE2 -nat_mort*fVE2 +m_mf*mVE2*(1-p_s) -m_fm_vax*fVE2
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1 +m_mf*mVI1*(1-p_s) -m_fm_vax*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2 +m_mf*mVI2*(1-p_s) -m_fm_vax*fVI2
    dmS = -(mbet1*mS*mI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet2*mS*mI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet1*mS*mVI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet2*mS*mVI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow +v_hat*mV -nat_mort*mS -m_mf*mS +m_fm*fS
    dmE1 = (mbet1*mS*mI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow +(mbet1*mS*mVI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -sig*mE1 -nat_mort*mE1 -m_mf*mE1 +m_fm*fE1
    dmE2 = (mbet2*mS*mI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow +(mbet2*mS*mVI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -sig*mE2 -nat_mort*mE2 -m_mf*mE2 +m_fm*fE2
    dmI1 = sig*mE1 -gamm*(1-p_1)*mI1 -mort*p_1*mI1 -nat_mort*mI1 -m_mf*mI1 +m_fm*fI1
    dmI2 = sig*mE2 -gamm*(1-p_2)*mI2 -mort*p_2*mI2 -nat_mort*mI2 -m_mf*mI2 +m_fm*fI2
    dmR = gamm*(1-p_1)*mI1 +gamm*(1-p_2)*mI2 +gamm*mVI1 +gamm*mVI2 -nat_mort*mR -m_mf*mR +m_fm*fR
    dmV = -(mbet1*mV*mI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet2*mV*mI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet1*mV*mVI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -(mbet2*mV*mVI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -v_hat*mV -nat_mort*mV -m_mf*mV +m_fm_vax*fV
    dmVE1 = (mbet1*mV*mI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow +(mbet1*mV*mVI1)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -sig*mVE1 -nat_mort*mVE1 -m_mf*mVE1 +m_fm_vax*fVE1
    dmVE2 = (mbet2*mV*mI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow +(mbet2*mV*mVI2)/(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)^N_pow -sig*mVE2 -nat_mort*mVE2 -m_mf*mVE2 +m_fm_vax*fVE2
    dmVI1 = sig*mVE1 -gamm*mVI1 -nat_mort*mVI1 -m_mf*mVI1 +m_fm_vax*fVI1
    dmVI2 = sig*mVE2 -gamm*mVI2 -nat_mort*mVI2 -m_mf*mVI2 +m_fm_vax*fVI2
    dM_FM = -M_FM +(1-p_s)*(m_fm*fS +m_fm*fE1 +m_fm*fE2 +m_fm*fI1 +m_fm*fI2 +m_fm*fR +m_fm_vax*fV +m_fm_vax*fVE1 +m_fm_vax*fVE2 +m_fm_vax*fVI1 +m_fm_vax*fVI2)
    dM_MF = -M_MF +m_mf*mS +m_mf*mE1 +m_mf*mE2 +m_mf*mI1 +m_mf*mI2 +m_mf*mR +m_mf*mV +m_mf*mVE1 +m_mf*mVE2 +m_mf*mVI1 +m_mf*mVI2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR, dfV, dfVE1, dfVE2, dfVI1, dfVI2,
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR, dmV, dmVE1, dmVE2, dmVI1, dmVI2,
                  dM_FM, dM_MF)))})}

# Set model 4 specific parameters and functions --------------------------------
c1 = 0.05
c2 = 0.4
N_pow = 1
#virs <- seq(1, 100, 1)
#betas <- c1*(virs)^c2
#plot(virs, betas, type='l')
# Initial susceptible population in farms
fS_init = (pop_size * 1/2) - 1
# Initial susceptible population in markets
mS_init = (pop_size * 1/2) - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Initial strain 1 infectious population in markets
mI1_init = 1
# Percent of susceptible chickens vaccinated in each time period
perc_vax = 0.26
# Time that perc_vax is vaccinated
inter_vax_time = 120 
# Vaccination rate of chickens of farms per susceptible chicken of farm per day
v = perc_vax / inter_vax_time 
# Rate of loss of immunity due to vaccination per chicken per day
v_hat = 1 / 120
# Percent sold in interval
perc_sold_per_farm = 0.1
# Days between successive sales of chickens of a farm
inter_sell_time_per_farm = 120
# Migration rate of chickens from farms to markets per chicken per day, if unvaccinated
m_fm = perc_sold_per_farm / inter_sell_time_per_farm
# Migration rate of chickens from markets to farms per chicken per day
m_mf = 1 / 7
# Ratio of contact rate in markets vs. farms
bet_mf_ratio = 5
# Threshold value for extinction
threshold_extinction = 2.2
# Percentage of market chickens that are to be immediately slaughtered
p_s = 0.8

# Assign model 4 specific equation and test_invade -----------------------------
eqn <- eqn_mod4_freq

# Run model --------------------------------------------------------------------
res_vir = 80
invade_vir = 80

# Strain specific virulence parameters
fbet1 <- (c1 * (res_vir)^c2)
mbet1 <- fbet1 * bet_mf_ratio
fbet2 <- (c1 * (invade_vir)^c2)
mbet2 <- fbet2 * bet_mf_ratio
p_1 <- ((res_vir) / 100)
p_2 <- ((invade_vir) / 100)

# Parameters
parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort, 
                v=v, v_hat=v_hat, m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf, mbet1=mbet1, mbet2=mbet2, p_s=p_s)

# Run resident strain until equilibrium
init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0, fV=0, fVE1=0, fVE2=0, fVI1=0, fVI2=0,
          mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=0, mR=0, mV=0, mVE1=0, mVE2=0, mVI1=0, mVI2=0,
          M_FM=0, M_MF=0)
time_eq1 <- seq(0, t_max_eq1, by = 1)
out_eq1 <- ode(y=init, times=time_eq1, eqn, parms=parameters)
out_eq1.df <- as.data.frame(out_eq1)

# morts <- seq(1, 100, 1)
# virs <- c1 * (morts) ^ c2
# plot(morts, virs, type='l')

# For debugging purposes
plot.out.df.mod4(out_eq1.df)

# Question 1: 
# Constant parameters
prev_f = 0.2
beta = 1
gamm = 1 / 7
N_pow = 1
N_f = 8e5
m_fm = 0.1 / 120
m_mf = 1 / 7

# Initial conditions
prev_m = 0.4
prop_recov = 0.3
N_m = (N_f * 0.01)
N_m = (1 - prop_recov) * N_m
recov_t = recov_0 = prop_recov * N_m

# Run loop
times = seq(1, 1000, 1)
prev_ms = c()
N_ms = c()
recov_ts = c()
for (t in times) {
  prev_ms = c(prev_ms, prev_m)
  N_ms = c(N_ms, N_m)
  N_m = N_m - m_mf*N_m + m_fm*N_f
  trans = (beta / (N_m)^N_pow)*((N_m)*(prev_m))*(N_m*(1-prev_m) - recov_t)
  prev_m_num = ((1 - m_mf) - gamm) * (N_m * prev_m) + m_fm*(N_f * prev_f) + trans
  prev_m_denom = ((1 - m_mf) + gamm) * (N_m * (1 - prev_m)) + m_fm * (N_f * (1 - prev_f)) - trans
  prev_m = prev_m_num / prev_m_denom
  recov_t = recov_t + gamm*(N_m * prev_m) - m_mf * (recov_t)
  recov_ts = c(recov_ts, recov_t)
}
plot(times, prev_ms, type='l')
plot(times, N_ms, type='l', ylim=c(0, 10000))
lines(times, recov_ts, type='l', col='red')

# Look at how the prevalence changes depending on each of the parameters?
# How many dimensions is it: (1) initial N_m (2) initial prop_recov (3) prev_f (4) m_fm (5) transmission (cleanliness)
# Really three dimensions if we get rid of initial conditions?
# Maybe just focus on the first two: m_fm and transmission which are control measures

# This model seems to show that you should find higher % recovered individuals in the markets than in farms, 
# by virtue of there not being as many susceptibles in markets compared to farms

# In this model so far, because there are no births or deaths in the system, if we solely look at markets, then there
# would be a greater benefit for more and more virulent strategies, since we do not feel any costs

# This is a frequency dependent model, looking at the relative importance of beta and
# prevalence of farms in determining the prevalence in markets
prev_fs <- seq(0.15, 0.4, 0.01)
betas <- seq(0, 5, 0.1)

res <- c()
recov_res <- c()
for (i in prev_fs) {
  for (j in betas) {
    print(paste0(i, '_', j))
    # Input parameters
    prev_f = prev_m = i
    beta <- j
    # Constant parameters
    gamm = 1 / 10
    N_pow = 0.2
    N_f = 1e6
    m_fm = 0.1 / 120
    m_mf = 1 / 7
    # Run loop with initial conditions
    prop_recov = 0.3
    N_m = (N_f * 0.01)
    curr_recov = prop_recov * N_m
    num_infect = prev_m * N_m
    num_noninfect = (1 - prev_m) * N_m
    mod_eqn <- function(time, state, parameters){
      with(as.list(c(state, parameters)),{
        dI = -(m_mf + gamm) * I + m_fm * (N_f * prev_f) + (beta / (I + NI)^N_pow) * (I) * (NI - R)
        dNI = -(m_mf) * NI + m_fm * (N_f * (1 - prev_f)) - (beta / (I + NI)^N_pow) * (I) * (NI - R)
        dR = (gamm * I) - (m_mf * R)
        return(list(c(dI, dNI, dR)))})}
    parameters <- c(m_mf = m_mf, gamm = gamm, m_fm = m_fm, beta = beta, N_f = N_f, prev_f = prev_f, N_pow = N_pow)
    init <- c(I = num_infect, NI = num_noninfect, R = curr_recov)
    times <- seq(1, 1000, 1)
    out <- ode(y=init, times=times, mod_eqn, parms=parameters)
    out.df <- as.data.frame(out)
    # plot(out.df$time, (out.df$I + out.df$NI), type='l', col='black', lwd=2, ylim=c(0, max(out.df$I + out.df$NI)))
    # lines(out.df$time, out.df$I, col='red')
    # lines(out.df$time, out.df$NI, col='green')
    # lines(out.df$time, out.df$R, col='blue')
    prev_m <- out.df$I[nrow(out.df)] / (out.df$I[nrow(out.df)] + out.df$NI[nrow(out.df)])
    res <- c(res, prev_m)
    recov <- out.df$R[nrow(out.df)] / (out.df$I[nrow(out.df)] + out.df$NI[nrow(out.df)])
    recov_res <- c(recov_res, recov)
  }
}
# I heat mat
res_matrix <- matrix(res, nrow=length(betas), ncol=length(prev_fs), byrow=F)
colnames(res_matrix) <- prev_fs
rownames(res_matrix) <- betas
melted <- melt(res_matrix)
colnames(melted) <- c('prev_f', 'beta', 'value')
temp_plot <- ggplot(melted, aes(y = prev_f, x = beta, fill = value)) + geom_tile() +
  xlab('prevalence of farms') + ylab('betas (per I per day hazard)') + ggtitle('prevalence in markets after 3 years') +
  theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
temp_plot

# R heat mat
res_matrix <- matrix(recov_res, nrow=length(betas), ncol=length(prev_fs), byrow=F)
colnames(res_matrix) <- prev_fs
rownames(res_matrix) <- betas
melted <- melt(res_matrix)
colnames(melted) <- c('prev_f', 'beta', 'value')
temp_plot <- ggplot(melted, aes(y = prev_f, x = beta, fill = value)) + geom_tile() +
  xlab('prevalence farms') + ylab('betas (per I per day hazard)') + ggtitle('Seroprevalence after 3 years') +
  theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
temp_plot

# Seems relevant that because there are no births in the market patch, there are always less susceptibles
# in the market patch compared to the farm patch

