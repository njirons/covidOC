rm(list=ls())

# format data and output of SEIRD model for input to NPI regression model

library(tidyr)
library(dplyr)
library(lubridate)
library(rjson)
library(readr)
library(cmdstanr)
options(mc.cores = parallel::detectCores())

# read in data
stan_data <- readRDS(file="data/stan_data.RData")
model_data <- readRDS(file="data/model_data.RData")
state_nums <- 1:stan_data$s

# randomly sample 100 posterior R0, E, I, and R_D trajectories from fitted models
n.traj <- 100 # number of posterior trajectories to use
n.samples <- 4000 # number of posterior samples
set.seed(123)
inds <- sample(1:n.samples,n.traj,replace=TRUE)
save(inds,file='data/traj_inds.Rda')
r0 <- list()
sus <- list()
exps <- list()
inf <- list()
rem_s <- list()
rem_d <- list()
d <- list()
for(state in state_nums){
    print(state)
  
    # load samples
    abbrev <- model_data$abbrev[state]
    load(paste0('state-sims/sims_',model_data$abbrev[state],'-',state,'.Rda'))
    
    # subset to random subsample
    sims <- sims[inds,]
    
    # extract data
    r0[[state]] <- as.matrix(sims %>% select(starts_with("r0[")))
    
    sus[[state]] <- as.matrix(sims %>% select(starts_with("sus[")))
    exps[[state]] <- as.matrix(sims %>% select(starts_with("exps[")))
    inf[[state]] <- as.matrix(sims %>% select(starts_with("inf[")))
    rem_s[[state]] <- as.matrix(sims %>% select(starts_with("rem_s[")))
    rem_d[[state]] <- as.matrix(sims %>% select(starts_with("rem_d[")))
    d[[state]] <- matrix(pmax(0,1 - 
      (
        sus[[state]] + exps[[state]] + inf[[state]] + rem_s[[state]] + rem_d[[state]]
      )
    ),nrow=n.traj)
}

# format data for NPI regression model
log_rt <- c()
sus_full <- c()
exps_full <- c()
inf_full <- c()
rem_s_full <- c()
rem_d_full <- c()
d_full <- c()
for(state in state_nums){
  print(state)
  
  log_rt <- cbind(log_rt, log(r0[[state]]))
  
  sus_full <- cbind(sus_full,sus[[state]])
  exps_full <- cbind(exps_full,exps[[state]])
  inf_full <- cbind(inf_full,inf[[state]])
  rem_s_full <- cbind(rem_s_full,rem_s[[state]])
  rem_d_full <- cbind(rem_d_full,rem_d[[state]])
  d_full <- cbind(d_full,d[[state]])
}

u_full <- matrix(unlist(stan_data$u),ncol = stan_data$p, byrow = TRUE)
p <- stan_data$p-1
u_full <- u_full[,1:p]
u <- matrix(NA,nrow=stan_data$n_error_terms,ncol=p)
t <- c(0,stan_data$error_term_inds[stan_data$t])
n_d <- diff(t)
exps_week <- c()
inf_week <- c()
rem_d_week <- c()
for(k in 1:stan_data$n_error_terms){
  week_inds <- which(stan_data$error_term_inds == k)
  if(length(week_inds) == 1){
    u[k,] <- u_full[week_inds,]  
  }else{
    u[k,] <- colMeans(u_full[week_inds,])  
  }  
  
  if(k %in% (t+1)){
    state <- which(t == k-1)
    exps_week <- cbind(exps_week,inf[[state]][,1])
    inf_week <- cbind(inf_week,rem_d[[state]][,1])
    rem_d_week <- cbind(rem_d_week,d[[state]][,1])
  }else{
    week_inds <- which(stan_data$error_term_inds == k-1)
    if(length(week_inds) == 1){
      exps_week <- cbind(exps_week,exps_full[,week_inds])
      inf_week <- cbind(inf_week,inf_full[,week_inds])
      rem_d_week <- cbind(rem_d_week,rem_d_full[,week_inds])
    }else{
      exps_week <- cbind(exps_week,rowSums(exps_full[,week_inds]))
      inf_week <- cbind(inf_week,rowSums(inf_full[,week_inds]))
      rem_d_week <- cbind(rem_d_week,rowSums(rem_d_full[,week_inds]))
    }
  }
}
state_inds <- c()
for(state in state_nums){
  state_inds <- c(state_inds,rep(which(state_nums == state),n_d[state]))
}
u_mean <- colMeans(u)

include_trunc <- 1
bhm_data <- list()
for(n in 1:n.traj){
  X <- cbind(
    stan_data$d_mult * rem_d_week[n,],
    stan_data$c_mult * inf_week[n,],
    stan_data$c_mult * exps_week[n,]
  )
  nX <- ncol(X)
  K <- p + nX + 1
  s <- stan_data$s
  
  bhm_data[[n]] <- list(
    p = p,
    nX = nX,
    K = K,
    b = rep(0,K),
    b_type = c(0,rep(-include_trunc,p),rep(0,nX)),
    s = s,
    t = t, 
    u = u, 
    u_mean = u_mean,
    X_all = X,
    X = X[,1:nX],
    log_rt = log_rt[n,],
    state_inds = state_inds,
    log_r0_upper = log(stan_data$r0_upper),
    r0_upper = stan_data$r0_upper,
    r0_lower = stan_data$r0_lower,
    r0_median = stan_data$r0_median,
    r0_mean = stan_data$r0_mean,
    r0_sd = stan_data$r0_sd
  )
}
saveRDS(bhm_data, file="data/bhm_data.RData")
