# marketVirEvol

Code to understand selection pressures in live poultry markets. Accompanies Sheen, Rasambainarivo, Saad-Roy, Grenfell, Metcalf, "Markets as drivers of selection for highly virulent poultry pathogens," Nature communications (2024)

## code
#### general: file containing functions used by files of code/main, as well as a testing file
- c1_c2_bounds.R: calculates range of c1 and c2 to test
- gen_fuctions.R: functions to calculate R0 when there is and is not differential migration of infectious poultry
- test_invade.R: testing file showing selection for lower virulence with higher cleaning

#### main: files used to create main results
- mod.R: calculate effects of cleaning and turnover rate on virulence evolution
- mod_all.R: calculate the cumulative effects of cleaning and turnover rate on virulence evolution
- mod_change.R: calculate the changes in virulence due to cleaning and turnover rate (Table S2)
- mod_all_change: calculate the cumulative changes in virulence due to cleaning and turnover rate (Table S2)
- ...the four scripts with "diff_migrate" in the title are the same as above, but for the sensitivity analysis of differential migration of infectious poultry (Table S3)

#### additional_figs: files to create Fig. 3 and Fig. S1
- fig3.R: script used to create Fig. 3 illustrating, for an example tradeoff curve, how persistence and turnover can select for higher virulece
- pip_fig/test_impulse2.R: script used to create pairwise invasibility results
- pip_fig/create_revision_fig.R: script used to plot pairwise invasibility plots
- pip_fig/high(low)Corrected.csv: pairwise invasibility result with high(low) cleaning, before rerunning for longer inconclusive cells
- pip_fig/kapHigh(Low)Mat2.csv: results of certain inconclusive cells of high(low) cleaing when running for longer
- pip_fig/kapHigh(Low)Mat2_corrected.csv: pairwise invasibility result with with high(low) cleaning, after rerunning for longer inconclusive cells
- pip_fig/out.png: pairwise invasibility results for the sensitivity analysis of our model when assuming distinct cohort periods (Fig. S1)

## code_output
#### obj: .RData objects created by the files in code/main
- psi(mm)_dens_shed_mod_X_Y_.RData: cleaning(turnover) result of runing mod.R for phi Y
- psi(mm)_dens_shed_diff_migrate_X_Y_.RData: cleaning(turnover) result of runing diff_migrate.R for phi Y
- dens_shed(_diff_migrate)_all_Y.RData: cumulative result from mod_all.R (diff_migrate_all.R) for phi Y

These programs are a work in progress, as we work to improve usability, error-catching, and speed of analysis. If you find errors, please contact Justin Sheen at jsheen (at) princeton (dot) edu.

Last Update: January 20, 2024
