# calibrating_large_scale_observations
Code and plots from simulations in the paper "Transporting treatment effects by calibrating large-scale observational outcomes"

The file 1d_sim.Rmd reproduces the univariate simulations (and corresponding figures in the paper), while the file multivar_sim.Rmd reproduces the multivariate simulations (and corresponding figures in the paper). Be sure to change the variable "PATH" at the top of each file to run the code on your machine.

Since the simulations take many hours to run, the code that actually runs the simulations has been commented out, the simulation output has been saved as various .rds files, and the code currently reads in these .rds files as the results of the sims to reproduce the plots in the paper.

