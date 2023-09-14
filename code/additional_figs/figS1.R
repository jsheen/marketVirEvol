

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
new_R0s_3 = c()
new_R0s_5 = c()
new_R0s_7 = c()
new_R0s_9 = c()
res = c()
env_orig = c()
env_3 = c()
env_5 = c()
env_7 = c()
env_9 = c()
for (alpha in alphas) {
  beta = c1 * (alpha ^ c2)
  lambda = beta * phi
  delta = alpha * c3
  
  get_R0 <- function(kappa) {
    temp = (epsilon * beta * lambda) / (kappa * psi * (m + gamma + delta + mu))
    temp
    orig_R0 = get_R0_shed(beta=beta, m_f=m_f, S_f=S_f, b=b, p=p, epsilon=epsilon, lambda=lambda,
                          sigma=sigma, nat_mort=mu, m=m, gamma=gamma, mort=delta, kappa=kappa, psi=psi, both=T)
    return(c(orig_R0, temp))
  }
  orig_R0s = c(orig_R0s, get_R0(kappa=1)[1])
  new_R0s_3 = c(new_R0s_3, get_R0(kappa=3)[1])
  new_R0s_5 = c(new_R0s_5, get_R0(kappa=5)[1])
  new_R0s_7 = c(new_R0s_7, get_R0(kappa=7)[1])
  new_R0s_9 = c(new_R0s_9, get_R0(kappa=9)[1])
  
  # This should be the difference between the two R0s i.e. reduction of second term due to cleaning
  temp = get_R0(kappa=1)[2]
  temp2 = get_R0(kappa=9)[2]
  to_add = ((temp - temp2) * (sigma / (sigma + mu + m)) * ((m_f * S_f) / (mu + m)))
  res = c(res, to_add) 
  
  # The following is just to get the environmental contribution of the original model
  get_env_res <- function(temp) {
    to_add2 = ((temp) * (sigma / (sigma + mu + m)) * ((m_f * S_f) / (mu + m)))
    return(to_add2)
  }
  env_orig = c(env_orig, get_env_res(get_R0(kappa=1)[2]))
  env_3 = c(env_3, get_env_res(get_R0(kappa=3)[2]))
  env_5 = c(env_5, get_env_res(get_R0(kappa=5)[2]))
  env_7 = c(env_7, get_env_res(get_R0(kappa=7)[2]))
  env_9 = c(env_9, get_env_res(get_R0(kappa=9)[2]))
}
plot(alphas, orig_R0s, type='l')
lines(alphas, new_R0s_9, type='l', col='red')
lines(alphas, new_R0s_9 + res, type='l') # just to double check
plot(alphas, res, type='l', col='blue')
plot(alphas, res / orig_R0s, type='l', col='purple')
all(res== cummax(res))

# Figure S1
## Panel A
par(mar=c(5,6,4,1)+.1)
plot(alphas, (env_orig / orig_R0s) * 100, type='l', col='black', xlab='α', ylab=expression('% of R'[0]*' is environment trans.'), main='(A)', lwd=3, cex.lab=1.5, ylim=c(0, max((env_orig / orig_R0s)*100)))
lines(alphas, (env_3 / new_R0s_3) * 100, type='l', col='purple', lwd=3)
lines(alphas, (env_5 / new_R0s_5) * 100, type='l', col='blue', lwd=3)
lines(alphas, (env_7 / new_R0s_7) * 100, type='l', col='red', lwd=3)
lines(alphas, (env_9 / new_R0s_9) * 100, type='l', col='orange', lwd=3)
legend(0, 12, legend=c("κ=1", "κ=3", 'κ=5', 'κ=7', 'κ=9'),
       col=c("black", "purple", 'blue', 'red', 'orange'), lty=1, cex=0.8, lwd=3)

## Panel B
par(mar=c(5,6,4,1)+.1)
plot(alphas, (env_orig / orig_R0s) * 100 - (env_9 / new_R0s_9) * 100, type='l', col='black', xlab='α', ylab=expression('% of R'[0]*' lost when κ=1 → κ=9'), main='(B)', lwd=3, lty='dotted',cex.lab=1.5, ylim=c(0, max((env_orig / orig_R0s)*100)))


# Change in ESS
which(orig_R0s == max(orig_R0s))
which(new_R0s_9 == max(new_R0s_9))


