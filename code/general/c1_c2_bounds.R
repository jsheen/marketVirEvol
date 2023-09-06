# Dens. dependent
virs <- seq(0, 1000, 1) # Should infect at least one per day
DFE_markets <- 2267.22
c2 <- 0.1
c1 <- 1 / 45
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect no more than 200 poultry initially
DFE_markets <- 2267.22
c2 <- 0.1
c1 <- 1 / 22
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect at least one per day
DFE_markets <- 2267.22
c2 <- 1
c1 <- 1 / 2300000
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)

virs <- seq(0, 1000, 1) # Should infect no more than 200 poultry initially
DFE_markets <- 2267.22
c2 <- 1
c1 <- 1 / 11300
infects <- c1 * (virs)^c2 * DFE_markets
plot(virs, infects, type='l')
max(infects)


# Freq. dependent
virs <- seq(0, 1000, 1)
c1 <- 50
c2 <- 0.1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 50
c2 <- 1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 0.51
c2 <- 0.1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')

c1 <- 0.1
c2 <- 1
toplot_betas <- c1 * (virs)^c2
plot(virs, toplot_betas, type='l')