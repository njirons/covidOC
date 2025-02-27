rm(list=ls())
library(rjson)

# read in python JSON data
stan_data <- fromJSON(file="data/stan_data.json")
model_data <- fromJSON(file="data/model_data.json")

# set epi parameters
stan_data$npi_cutoffs <- rep(0,stan_data$s)
stan_data$gamma_mean <- 5
stan_data$delta_mean <- 5.5
stan_data$mu_mean <- 10.5
stan_data$init_prop <- 0.05
stan_data$init_props <- rep(stan_data$init_prop,stan_data$s)
stan_data$r0_upper <- 6.5
stan_data$kappa_d_upper <- 1000
stan_data$ccd_sd <- 4.116
stan_data$ccd_lower <- stan_data$delta_mean
stan_data$ccd_upper <- stan_data$delta_mean + 
  stan_data$gamma_mean + stan_data$mu_mean 
stan_data$ccd_mean <- stan_data$delta_mean +
  stan_data$gamma_mean + stan_data$mu_mean - 8.053
stan_data$ifr_lower <- stan_data$ifr_mean - 4*stan_data$ifr_sd
stan_data$ifr_upper <- stan_data$ifr_mean + 4*stan_data$ifr_sd
stan_data$n_compartments <- 5

stan_data$d_mult <- 1000
stan_data$c_mult <- 10

# save final data set
saveRDS(stan_data, file="data/stan_data.RData")
saveRDS(model_data, file="data/model_data.RData")
