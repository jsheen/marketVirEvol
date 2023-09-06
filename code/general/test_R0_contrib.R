

source('~/marketVirEvol/code/general/gen_functions.R')

b = 0
p = 0
m_f = 0.1 / 120
N_f = 1e6
prev_f = 0.12
I_f = E_f = N_f * (prev_f / 2)
seroprev_f <- 0.402
R_f = seroprev_f * (E_f + (1 - prev_f) * N_f)
S_f = ((1 - prev_f) * N_f) - R_f
c1 = 1 / 2300
c2 = 0.6
c3 = 1 / 1000
epsilon = 1
phi = 1
psi = 1 / 5
m = 1 / 5.5
gamma = 1 / 5
mu = 1 / 365
sigma = 1 / 5

alphas = seq(1, 1000, 1)
orig_R0s = c()
new_R0s = c()
res = c()
res2 = c()
for (alpha in alphas) {
  beta = c1 * (alpha ^ c2)
  lambda = beta * phi
  delta = alpha * c3
  
  kappa = 1
  temp = (epsilon * beta * lambda) / (kappa * psi * (m + gamma + delta + mu))
  temp
  orig_R0 = get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                   sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=delta, kappa=kappa, psi=psi, both=T)
  orig_R0s = c(orig_R0s, orig_R0)
  
  kappa = 9
  temp2 = (epsilon * beta * lambda) / (kappa * psi * (m + gamma + delta + mu))
  temp2
  new_R0 = get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                        sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=delta, kappa=kappa, psi=psi, both=T)
  new_R0s = c(new_R0s, new_R0)
  
  to_add = ((temp - temp2) * (sigma / (sigma + mu + m)) * ((m_f * S_f) / (mu + m)))
  res = c(res, to_add) # This should be the difference between the two R0s i.e. reduction of second term due to cleaning
  
  # The following is just to get the environmental contribution of the original model
  to_add2 = ((temp) * (sigma / (sigma + mu + m)) * ((m_f * S_f) / (mu + m)))
  res2 = c(res2, to_add2)
}
plot(alphas, orig_R0s, type='l')
lines(alphas, new_R0s, type='l', col='red')
lines(alphas, new_R0s + res, type='l') # just to double check
plot(alphas, res, type='l', col='blue')
plot(alphas, res / orig_R0s, type='l', col='purple')
all(res== cummax(res))

# Figure S2
plot(alphas, (res2 / orig_R0s) * 100, type='l', col='black', xlab='α', ylab='% of R0 is environment trans.', lwd=5, cex.lab=1.5)

# Change in ESS
which(orig_R0s == max(orig_R0s))
which(new_R0s == max(new_R0s))


