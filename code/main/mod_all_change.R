# Script in order to get the percent of tradeoff curves with zero change in global ESS
# and percent of tradeoff curves with nonzero change in global ESS, as well as the maximum
# change in ESS, all due to either increases in turnover rate, or decreases in cleaning.

# When there is a relationship between lambda and virulence --------------------
get_perc_nonzero_change_all <- function(phi) {
  load(paste0("~/marketVirEvol/code_output/obj/dens_shed_all_", phi, ".RData"))
  max_diffs <- c()
  zero_diffs <- 0
  nonzero_diffs <- 0
  for (diff_vir in diff_virs) {
    # All differences should be positive, return error if not
    if (sum(diff_vir < 0) != 0) {
      stop('Error, no values should be negative.')
    }
    # Check that it monotonically increases from left to right
    for (i in 1:nrow(diff_vir)) {
      foc_vec <- diff_vir[i,]
      if(!all(foc_vec == cummax(foc_vec))){
        stop('Error, not monotonically increasing from left to right.')
      }
    }
    # Check that it monotonically decreases from top to bottom
    for (j in 1:ncol(diff_vir)) {
      foc_vec <- diff_vir[,j]
      if(!all(foc_vec == cummin(foc_vec))){
        stop('Error, not monotonically decreasing from top to bottom.')
      }
    }
    
    max_diffs <- c(max_diffs, max(diff_vir))
    if (length(which(diff_vir > 0)) > 0) {
      nonzero_diffs <- nonzero_diffs + 1
    } else {
      zero_diffs <- zero_diffs + 1
    }
  }
  print(paste0('N: ', length(diff_virs)))
  print(paste0('% nonzero_diffs: ', (nonzero_diffs / (zero_diffs + nonzero_diffs)) * 100))
  hist(max_diffs) # Max diff for each tradeoff curve
  max(max_diffs)
}
get_perc_nonzero_change_all(phi=0.1)
get_perc_nonzero_change_all(phi=1)
get_perc_nonzero_change_all(phi=10)
get_perc_nonzero_change_all(phi=100)


