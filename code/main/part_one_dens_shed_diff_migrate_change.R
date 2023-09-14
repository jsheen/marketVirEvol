# When there is a relationship between lambda and virulence and migrate rate differs -------
lambda = 0.0434782608695652
phi = 1
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_diff_migrate", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
constrain_R0 <- which(diff_virs_R0 <= 100 & diff_virs_R0 >= 1)
mm <- mm[constrain_R0]
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_diff_migrate", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
diff_virs_R0 <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
constrain_R0 <- which(diff_virs_R0 <= 100 & diff_virs_R0 >= 1)
psi <- psi[constrain_R0]
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
# Q0 result: true, there is a single optimum in this model for parameters tested
all(opt_psi_res)
# Q1 result: true, as psi increases, R0 decreases for all virulence strategies
all(R0_psi_res)
# Q2 result: false, though it is generally true that as psi increases, virulence strategies decrease a lot after optimum. This makes sense, as as psi increases, the very virulent strategies suffer the most by not having as great transmission gains.
all(flat_psi_res)
# Q3 result: true, as psi increases, optimal virulence strategies decrease. This is because there is less transmission benefits.
all(inc_psi_res)
rm(list=ls())