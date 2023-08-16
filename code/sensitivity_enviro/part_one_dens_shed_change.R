
# 1)
lambda = 0.1
phi = 'NA'
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
rm(list=ls())

# 2)
lambda = 1
phi = 'NA'
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
rm(list=ls())

# 3)
lambda = 0.0434782608695652
phi = 1
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
rm(list=ls())

# 4)
lambda = 0.434782608695652
phi = 10
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
rm(list=ls())

# 5)
lambda = 4.34782608695652
phi = 100
load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_", lambda, "_", phi, ".RData"))
mm <- sapply(diff_virs_1, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(mm == 0)) / length(mm)
length(which(mm > 0)) / length(mm)
if (length(which(mm < 0)) > 0) {
  stop('Error.')
}
load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_", lambda, "_", phi, ".RData"))
psi <- sapply(diff_virs_2, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
length(which(psi == 0)) / length(psi)
length(which(psi > 0)) / length(psi)
if (length(which(psi < 0)) > 0) {
  stop('Error.')
}
rm(list=ls())


