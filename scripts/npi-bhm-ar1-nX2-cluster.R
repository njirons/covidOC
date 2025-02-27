rm(list=ls())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
traj.ind <- as.integer(args[1])

library(cmdstanr)
options(mc.cores = parallel::detectCores())

# read in data
stan_data <- readRDS(file="data/stan_data.RData")
model_data <- readRDS(file="data/model_data.RData")

# format input data
file.copy('scripts/npi-bhm-ar1-robust.stan',paste0('scripts/npi-bhm-ar1-robust-',traj.ind,'.stan'))
file <- file.path(paste0('scripts/npi-bhm-ar1-robust-',traj.ind,'.stan')) 
mod <- cmdstan_model(file)
file.remove(file)
num.chains <- 4

bhm_data <- readRDS('data/bhm_data.RData')
nX <- 2
include_trunc <- 1 # truncate NPI prior

bhm_data[[traj.ind]]$nX <- nX
bhm_data[[traj.ind]]$X <- matrix(bhm_data[[traj.ind]]$X_all[,1:nX],ncol=nX)
p <- bhm_data[[traj.ind]]$p
bhm_data[[traj.ind]]$K <- p + nX + 1
s <- bhm_data[[traj.ind]]$s
K <- bhm_data[[traj.ind]]$K
n_ar <- bhm_data[[traj.ind]]$n_ar
bhm_data[[traj.ind]]$b <- rep(0,K)
bhm_data[[traj.ind]]$b_type <- c(0,rep(-include_trunc,p),rep(0,nX))

init <- list()
for(j in 1:num.chains){
  init[[j]] <- list(
    u_s = matrix(0.5,nrow=s,ncol=K),
    phi_ar = rep(0,n_ar),
    sigma = 0.4,
    sd_re = rep(0.4,K),
    L_corr = diag(K),
    mu = c(log(3),rep(-0.1,p),rep(-0.3,nX))
  )  
}

fit <- mod$sample(
  data = bhm_data[[traj.ind]],
  seed = 123,
  chains = num.chains, parallel_chains = num.chains,
  iter_warmup = 2000, iter_sampling = 25,
  max_treedepth = 10, adapt_delta=0.8, save_warmup=FALSE,
  init = init
) 

sims <- fit$draws(format='df')

saveRDS(fit, file=paste0('npi-bhm-fits/npi-bhm-ar1-nX',nX,'-',traj.ind,'.rds'))
save(sims, file=paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-',traj.ind,'.Rda'))
