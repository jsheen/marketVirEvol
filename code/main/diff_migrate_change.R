# Script in order to get the percent of tradeoff curves with zero change in global ESS
# and percent of tradeoff curves with nonzero change in global ESS, as well as the maximum
# change in ESS, all due to either increases in turnover rate, or decreases in cleaning.

# When there is a relationship between lambda and virulence --------------------
get_perc_nonzero_change <- function(lambda, phi) {
  max_R0_consider = 100
  
  load(paste0("~/marketVirEvol/code_output/obj/mm_dens_shed_diff_migrate_", lambda, "_", phi, ".RData"))
  mm <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
  diff_virs_R0_mm <- sapply(diff_virs_m_m, function(x) as.numeric(strsplit(x, ',')[[1]][5]))
  constrain_R0_mm <- which(diff_virs_R0_mm <= max_R0_consider & diff_virs_R0_mm >= 1)
  mm <- mm[constrain_R0_mm]
  print(paste0('N (mm): ', length(mm)))
  print(paste0('prop. zero change (mm): ', length(which(mm == 0)) / length(mm)))
  print(paste0('prop. nonzero change (mm): ',  length(which(mm > 0)) / length(mm)))
  print(paste0('max nonzero change (mm): ', max(mm)))
  if (length(which(mm < 0)) > 0) {
    stop('Error.')
  }
  load(paste0("~/marketVirEvol/code_output/obj/psi_dens_shed_diff_migrate_", lambda, "_", phi, ".RData"))
  psi <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][4]))
  diff_virs_R02_psi <- sapply(diff_virs_psi, function(x) as.numeric(strsplit(x, ',')[[1]][6]))
  constrain_R0_psi <- which(diff_virs_R02_psi <= max_R0_consider & diff_virs_R02_psi >= 1)
  psi <- psi[constrain_R0_psi]
  print(paste0('N (psi): ', length(psi)))
  print(paste0('prop. zero change (psi): ', length(which(psi == 0)) / length(psi)))
  print(paste0('prop. nonzero change (psi): ',  length(which(psi > 0)) / length(psi)))
  print(paste0('max nonzero change (psi): ', max(psi)))
  if (length(which(psi < 0)) > 0) {
    stop('Error.')
  }
  # Q0 result: true, there is a single optimum in this model for parameters tested
  if (!all(opt_psi_res)) {
    stop('Error (1).')
  }
  # Q1 result: true, as psi increases, R0 decreases for all virulence strategies
  if (!all(R0_psi_res)) {
    stop('Error (2).')
  }
  # Q2 result: false, though it is generally true that as psi increases, virulence strategies decrease a lot after optimum. This makes sense, as as psi increases, the very virulent strategies suffer the most by not having as great transmission gains.
  all(flat_psi_res)
  # Q3 result: true, as psi increases, optimal virulence strategies decrease. This is because there is less transmission benefits.
  if (!all(inc_psi_res)) {
    stop('Error (3).')
  }
}
get_perc_nonzero_change(lambda=0.00434782608695652, phi=0.1)
get_perc_nonzero_change(lambda=0.0434782608695652, phi=1)
get_perc_nonzero_change(lambda=0.434782608695652, phi=10)
get_perc_nonzero_change(lambda=4.34782608695652, phi=100)


