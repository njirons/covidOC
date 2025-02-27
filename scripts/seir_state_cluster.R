rm(list=ls())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
state_num <- as.integer(args[1])

library(rjson)
library(cmdstanr)
options(mc.cores = parallel::detectCores())

# read in data
stan_data <- readRDS(file="data/stan_data.RData")
model_data <- readRDS(file="data/model_data.RData")

# format input data
state_data <- list(
  n_compartments = stan_data$n_compartments,
  n_d = stan_data$n_d[state_num],
  d = stan_data$d[(stan_data$t_d[state_num]+1):stan_data$t_d[state_num+1]],
  c = stan_data$c[(stan_data$t_d[state_num]+1):stan_data$t_d[state_num+1]],
  d0 = stan_data$d0[state_num],
  c0 = stan_data$c0[state_num],
  n = stan_data$n[state_num],
  init_props = stan_data$init_prop,
  r0_upper = stan_data$r0_upper,
  r0_lower = stan_data$r0_lower,
  r0_mean = stan_data$r0_mean,
  r0_median = stan_data$r0_median,    
  r0_sd = stan_data$r0_sd,
  ifr_upper = stan_data$ifr_upper[state_num],
  ifr_lower = stan_data$ifr_lower[state_num],
  ifr_mean = stan_data$ifr_mean[state_num],
  ifr = stan_data$ifr_mean[state_num],
  ifr_sd = stan_data$ifr_sd[state_num],
  gamma = 1/stan_data$gamma_mean,
  delta = 1/stan_data$delta_mean,
  mu = 1/stan_data$mu_mean,
  n0_d = stan_data$n0_d[state_num],
  n0_c = stan_data$n0_c[state_num],
  nz_d = stan_data$nz_d[state_num],
  nz_c = stan_data$nz_c[state_num],
  d_nz_inds = 
    stan_data$d_nz_inds[
      (stan_data$nz_d_cumsum[state_num]+1):
        stan_data$nz_d_cumsum[state_num+1]] - 
    stan_data$t_d[state_num],
  c_nz_inds = 
    stan_data$c_nz_inds[
      (stan_data$nz_c_cumsum[state_num]+1):
        stan_data$nz_c_cumsum[state_num+1]] - 
    stan_data$t_d[state_num],
  d_z_inds = 
    stan_data$d_z_inds[
      (stan_data$n0_d_cumsum[state_num]+1):
        stan_data$n0_d_cumsum[state_num+1]] -
    stan_data$t_d[state_num],
  c_z_inds =
    stan_data$c_z_inds[
      (stan_data$n0_c_cumsum[state_num]+1):
        stan_data$n0_c_cumsum[state_num+1]] -
    stan_data$t_d[state_num],
  error_term_inds = stan_data$error_term_inds[
    (stan_data$t_d[state_num]+1):stan_data$t_d[state_num+1]],
  sigma_upper = stan_data$sigma_upper,
  alpha_sigma = stan_data$alpha_sigma,
  beta_sigma = stan_data$beta_sigma,
  kappa_d_upper = stan_data$kappa_d_upper,
  ccd_lower = stan_data$ccd_lower,
  ccd_upper = stan_data$ccd_upper,
  ccd_mean = stan_data$ccd_mean,
  ccd_sd = stan_data$ccd_sd
)
state_data$n0_c_inds <- length(state_data$c_z_inds)
state_data$n0_d_inds <- length(state_data$d_z_inds)
state_data$error_term_inds <- 
  state_data$error_term_inds - 
  min(state_data$error_term_inds) + 1
state_data$n_error_terms <- length(unique(state_data$error_term_inds))

save(state_data,file=paste0('state-data/data_',model_data$abbrev[state_num],'-',state_num,'.Rda'))

# build stan model
if(state_num %in% c(1,8,9,16,36,40,42,47,51)){# poisson dispersed data
  file.copy('scripts/seir-state-poisson.stan',paste0('scripts/seir-state-poisson',state_num,'.stan'))
  file <- file.path(paste0('scripts/seir-state-poisson',state_num,'.stan'))    
}else{# over-dispersed NB data
  file.copy('scripts/seir-state-nb.stan',paste0('scripts/seir-state-nb',state_num,'.stan'))
  file <- file.path(paste0('scripts/seir-state-nb',state_num,'.stan'))      
}
mod <- cmdstan_model(file)
file.remove(file)

# fit model
warmup <- 5000
iter <- 1000
chains <- 4
max_treedepth <- 10
adapt_delta <- 0.8
fit <- mod$sample(
  data = state_data,
  seed = 123,
  chains = chains, parallel_chains = chains,
  iter_warmup = warmup, iter_sampling = iter,
  max_treedepth = max_treedepth, adapt_delta = adapt_delta, 
  save_warmup=FALSE
)

# save model fit
fit$save_object(file = paste0('state-fits/fit_',model_data$abbrev[state_num],'-',state_num,'.rds'))

sims <- fit$draws(format = "df")
save(sims,file=paste0('state-sims/sims_',model_data$abbrev[state_num],'-',state_num,'.Rda'))
