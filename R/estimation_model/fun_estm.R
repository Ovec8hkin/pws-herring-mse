# fun_alt_assess.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Simplified assessment model that takes true numbers from operating model to produce assessment estimates

setwd("/Volumes/G-DRIVE mobile SSD R-Series/IDrive-Sync/Thesis/BASA_handoff/Thesis/ASA/2021/admb/base/")

fun_alt_assess <- function(pop_dyn, prev_assess_dir){
  
  library(tidyverse)
  
  list2env(pop_dyn)
  
  #### Read in parameter posteriors from previous full assessment
  pars <- read_csv("mcmc_out/iterations.csv",col_names = T)
  
  #### Sample parameter vectors from most recent MCMC - should we retain these same parameter vectors over all years in MSE?
  n_samp <- 100
  samp_ind <- sample(1:nrow(pars),n_samp)
  pars_samp <- pars[samp_ind,]
  
  #### Create vector of possible recruitment devs that spans historical range (length 100, or 1000)
  # Includes all posterior samples as possible values
  n_posrec <- 1000
  future_rec <- pars %>% select(starts_with("annual_age0devs")) %>%
    pivot_longer(everything()) %>%
    # summarise(annual_age0devs_next=runif(1000,min(value),max(value)))
    summarise(annual_age0devs_next=seq(min(value),max(value),length.out=n_posrec))
  
  #### Expand all combinations of parameter vectors AND unique possible recruitment values
  ind <- expand_grid(x=1:n_samp,y=1:n_posrec)
  pars_expanded <- bind_cols(pars_samp[ind$x,],future_rec[ind$y,])
  
  #### Use these "new" parameter vectors to calculate SSB and NAA (KEEP THESE), and predictions of data for fitting
  # How? 1) Take many of the equations from BASA in separate function to project out SSB over all years
  
  # Calculate new total posterior density based on likelihood functions used to determine fit for new data
  
  # Sum posterior densities of all vectors, and divide each posterior density by this sum to obtain posterior probabilities
  
  # Resample from all these vector using the posterior probabilities as weights, to generate new posterior distributions of SSB and NAA
  
  
  # OUTPUT the posterior of the estimated numbers-at-age 

  
  return(posterior_nya)
}