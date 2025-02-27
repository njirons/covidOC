rm(list=ls())

# Example to run policy optimization for a single (US state, cost function, NPI model) specification

state.num <- 1 # Alaska
cost.setting <- 'nuL-sL-wM' # baseline (low infection (nu) cost, low school (s) cost, medium work (w) cost)
nX <- 2 # model with nX confounders

library(cmdstanr)
library(rjson)
library(readr)
library(optimParallel)
library(snow)
library(dplyr)
options(mc.cores = parallel::detectCores())

source('scripts/oc-functions.R')

# load state data
stan.data <- readRDS(file="data/stan_data.RData")
model.data <- readRDS(file="data/model_data.RData")
cost.data <- read_csv('data/state-costs.csv')
df <- readRDS('data/bhm_data.RData')[[1]]

# extract relevant parameter posteriors
state.ind <- which(df$state_inds == state.num)
weeks <- state.ind
n.weeks <- length(weeks)
n.days <- stan.data$n_d[state.num]
week.ind <- 
  stan.data$error_term_inds[stan.data$state_inds==state.num] - min(weeks) + 1
n.pop <- stan.data$n[state.num]

## load SEIRD model samples
abbrev <- model.data$abbrev[state.num]
load(paste0('state-sims/sims_',abbrev,'-',state.num,'.Rda'))
sims.seird <- sims; rm(sims)

# load full NPI regression model samples
load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-full.Rda'))

## determine total posterior sample size as minimum of SEIRD and NPI models
n.iter <- min(dim(sims.seird)[1],dim(sims.npi)[1])

# subset samples to have same size
set.seed(123)
inds.seird <- sample(1:dim(sims.seird)[1],n.iter)
inds.npi <- sample(1:dim(sims.npi)[1],n.iter)
sims.seird.sub <- sims.seird[inds.seird,]
sims.npi.sub <- sims.npi[inds.npi,]
rm(sims.npi)

mu <- as.matrix(sims.npi.sub %>% select(starts_with('mu[')))
mu_s <- as.matrix(sims.npi.sub %>% select(starts_with(paste0('mu_s[',state.num,','))))
eps <- as.matrix(sims.npi.sub %>% select(starts_with('eps[')))[,weeks]
d0 <- sims.seird.sub$ifr * stan.data$init_props[state.num] *
  (sims.seird.sub %>% select(starts_with(paste0('x0[',stan.data$n_compartments))))[[1]]

# set cost parameters
cost.ind <- which(cost.data$state == abbrev)

cost.nu.high = mean(sims.seird.sub$ifr) * cost.data$vsl.high[cost.ind] +
  cost.data$cost.fear[cost.ind] + cost.data$cost.productivity[cost.ind] + 
  cost.data$cost.medical[cost.ind]
cost.nu.mid = mean(sims.seird.sub$ifr) * cost.data$vsl.mid[cost.ind] +
  cost.data$cost.fear[cost.ind] + cost.data$cost.productivity[cost.ind] + 
  cost.data$cost.medical[cost.ind]
cost.nu.low = mean(sims.seird.sub$ifr) * cost.data$vsl.low[cost.ind] +
  cost.data$cost.fear[cost.ind] + cost.data$cost.productivity[cost.ind] + 
  cost.data$cost.medical[cost.ind]

cost.work.low <- cost.data$cost.work[cost.ind]/2
cost.work.mid <- cost.data$cost.work[cost.ind]
cost.work.high <- cost.data$cost.work[cost.ind]*1.5

cost.social.low <- cost.data$cost.social[cost.ind]/2
cost.social.mid <- cost.data$cost.social[cost.ind]
cost.social.high <- cost.data$cost.social[cost.ind]*1.5

# run OC for various parameter settings
path <- paste0('oc-',cost.setting,'-sdM-nX',nX)
cost.social <- cost.social.mid
if(cost.setting == 'nuL-sL-wM'){
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.mid 
}else if(cost.setting == 'nuL-sL-wL'){
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.low   
}else if(cost.setting == 'nuL-sL-wH'){
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.high   
}else if(cost.setting == 'nuL-sH-wL'){
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.high[cost.ind]
  cost.work <- cost.work.low    
}else if(cost.setting == 'nuL-sH-wH'){
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.high[cost.ind]
  cost.work <- cost.work.high    
}else if(cost.setting == 'nuH-sL-wL'){
  cost.nu <- cost.nu.high
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.low    
}else if(cost.setting == 'nuH-sL-wH'){
  cost.nu <- cost.nu.high
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.high    
}else if(cost.setting == 'nuH-sH-wL'){
  cost.nu <- cost.nu.high
  cost.learning <- cost.data$cost.learning.high[cost.ind]
  cost.work <- cost.work.low    
}else if(cost.setting == 'nuH-sH-wH'){
  cost.nu <- cost.nu.high
  cost.learning <- cost.data$cost.learning.high[cost.ind]
  cost.work <- cost.work.high    
}else{
  stop("Not a valid cost setting.", call.=FALSE)
}

normal.school.year <- c(
  1, # winter break, first week of january
  rep(0,11),
  1, # spring break
  rep(0,11),
  rep(1,12), # summer break, mid june thru mid september
  rep(0,10),
  1, # thanksgiving break
  rep(0,4),
  1 # winter break, last week of december
)
normal.school.year <- tail(normal.school.year,n.weeks)

params <- list(
  n.iter = n.iter,
  nX = nX,
  mu0 = mu_s[,1],
  mu_npi = mu_s[,2:(df$p+1)],
  mu_X = mu_s[,(df$p+2):ncol(mu_s)],
  phi.ar = sims.npi.sub$phi_ar,
  sigma = sims.npi.sub$sigma,
  eps = eps,
  n.days = n.days,
  n.weeks = n.weeks,
  week.ind = week.ind,
  sus0 = sims.seird.sub$`sus[1]`,
  exps0 = sims.seird.sub$`exps[1]`,
  inf0 = sims.seird.sub$`inf[1]`,
  rem_s0 = sims.seird.sub$`rem_s[1]`,
  rem_d0 = sims.seird.sub$`rem_d[1]`,
  d0 = d0,
  ifr = sims.seird.sub$ifr,
  gamma = 1/stan.data$gamma_mean,
  delta = 1/stan.data$delta_mean,
  mu = 1/stan.data$mu_mean,
  n.pop = n.pop,
  u = t(df$u[state.ind,1:df$p]),
  u_mean = df$u_mean[1:df$p],
  d_mult = stan.data$d_mult,
  c_mult = stan.data$c_mult,
  cost.nu = cost.nu,
  cost.info = cost.data$cost.info[cost.ind],
  cost.mask = cost.data$cost.mask[cost.ind],
  cost.test = cost.data$cost.test[cost.ind],
  cost.trace = cost.data$cost.trace[cost.ind],
  cost.work = cost.work,
  cost.social = cost.social,
  cost.learning = cost.learning,
  cost.schools.direct = cost.data$cost.schools.direct[cost.ind],
  time.learning = cost.data$time.learning[cost.ind],
  weekly.discount = cost.data$weekly.discount[cost.ind],
  daily.discount = cost.data$daily.discount[cost.ind],
  normal.school.year = normal.school.year,
  n.school.weeks = sum(1-normal.school.year)
)

# define initializations
u.init <- rbind(
  runif(params$n.school.weeks+5*params$n.weeks),
  runif(params$n.school.weeks+5*params$n.weeks),
  runif(params$n.school.weeks+5*params$n.weeks),
  runif(params$n.school.weeks+5*params$n.weeks),
  rep(1,params$n.school.weeks+5*params$n.weeks),
  rep(0,params$n.school.weeks+5*params$n.weeks),
  rep(0.5,params$n.school.weeks+5*params$n.weeks),
  c(rep(0.1,params$n.school.weeks),rep(0.9,5*params$n.weeks)),  
)

u.oc <- matrix(NA,nrow=nrow(u.init),ncol=ncol(u.init))
val.oc <- rep(NA,nrow(u.init))

# find optimal control
for(j in 1:nrow(u.init)){
  oc.eval <- optimParallel(
    par = u.init[j,],
    fn = cost.final,
    method = 'L-BFGS-B',
    params = params,
    lower = 0, upper = 1,
    control = list(trace = 6),
    parallel = list(forward = TRUE)
  )
  save(oc.eval,file=paste0('oc-results/',abbrev,'-',path,'-',j,'.Rda'))
  
  u.oc[j,] <- oc.eval$par
  val.oc[j] <- oc.eval$value
}

save(u.oc,file=paste0('oc-results/',abbrev,'-',path,'-mat.Rda'))
write.csv(u.oc,file=paste0('oc-results/',abbrev,'-',path,'-mat.csv'))
save(val.oc,file=paste0('oc-results/',abbrev,'-',path,'-val.Rda'))
write.csv(val.oc,file=paste0('oc-results/',abbrev,'-',path,'-val.csv'))

val.opt <- which.min(val.oc)
u.oc.vec <- u.oc[val.opt,]
u.oc.school <- params$normal.school.year
u.oc.school[u.oc.school == 0] <- u.oc.vec[1:params$n.school.weeks]
u.oc.vec <- c(u.oc.school,u.oc.vec[-(1:params$n.school.weeks)])
u.oc.opt <- matrix(u.oc.vec,ncol=params$n.weeks,byrow=TRUE)
u.oc.opt <- make.npis(u.oc.opt,params)
save(u.oc.opt,file=paste0('oc-results/',abbrev,'-',path,'-opt.Rda'))
write.csv(u.oc.opt,file=paste0('oc-results/',abbrev,'-',path,'-opt.csv'))


