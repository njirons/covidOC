rm(list=ls())

# sensitivity analysis of optimal control policies

library(cmdstanr)
library(rjson)
library(readr)
library(optimParallel)
library(snow)
library(dplyr)
library(lubridate)
library(latex2exp)
library(vioplot)
options(mc.cores = parallel::detectCores())

source('scripts/oc-functions.R')

# load state data
stan.data <- readRDS(file="data/stan_data.RData")
model.data <- readRDS(file="data/model_data.RData")
cost.data <- read_csv('data/state-costs.csv')
df <- readRDS('data/bhm_data.RData')[[1]]

nX <- 3
load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-full.Rda'))
sims.npi <- sims; rm(sims)

paths <- c(
  'oc-nuH-sL-wH-sdM',
  'oc-nuH-sH-wH-sdM',
  'oc-nuH-sL-wL-sdM',
  'oc-nuH-sH-wL-sdM',
  'oc-nuL-sL-wH-sdM',
  'oc-nuL-sH-wH-sdM',
  'oc-nuL-sL-wL-sdM',
  'oc-nuL-sH-wL-sdM'
)

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

# to plot costs under various policies

# generate results

cost.oc <- matrix(NA,nrow=length(paths),ncol=51)
cost.open <- matrix(NA,nrow=length(paths),ncol=51)
cost.obs <- matrix(NA,nrow=length(paths),ncol=51)
cost.full <- matrix(NA,nrow=length(paths),ncol=51)
cost.noschool.full <- matrix(NA,nrow=length(paths),ncol=51)
cost.noschool.obs <- matrix(NA,nrow=length(paths),ncol=51)

u.oc <- list()

for(state.num in 1:51){
  print(state.num)
  
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
  
  ## determine total posterior sample size as minimum of SEIRD and NPI models
  n.iter <- min(dim(sims.seird)[1],dim(sims.npi)[1])
  
  # subset samples to have same size
  set.seed(123)
  inds.seird <- sample(1:dim(sims.seird)[1],n.iter)
  inds.npi <- sample(1:dim(sims.npi)[1],n.iter)
  sims.seird.sub <- sims.seird[inds.seird,]
  sims.npi.sub <- sims.npi[inds.npi,]
  
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
  cost.work.high <- cost.data$cost.work[cost.ind]*2
  
  cost.social.low <- cost.data$cost.social[cost.ind]/2
  cost.social.mid <- cost.data$cost.social[cost.ind]
  cost.social.high <- cost.data$cost.social[cost.ind]*2
  
  u.oc.tmp <- list()
  for(j in 1:length(paths)){
    
    # baseline cost function
    path <- paths[j]
    file.path <- paste0(path,'-nX',nX)
    if(path == 'oc-nuH-sL-wH-sdM'){
      cost.nu <- cost.nu.high
      cost.learning <- cost.data$cost.learning.low[cost.ind]
      cost.work <- cost.work.high
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuH-sH-wH-sdM'){
      cost.nu <- cost.nu.high
      cost.learning <- cost.data$cost.learning.high[cost.ind]
      cost.work <- cost.work.high
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuH-sL-wL-sdM'){
      cost.nu <- cost.nu.high
      cost.learning <- cost.data$cost.learning.low[cost.ind]
      cost.work <- cost.work.low
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuH-sH-wL-sdM'){
      cost.nu <- cost.nu.high
      cost.learning <- cost.data$cost.learning.high[cost.ind]
      cost.work <- cost.work.low
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuL-sL-wH-sdM'){
      cost.nu <- cost.nu.low
      cost.learning <- cost.data$cost.learning.low[cost.ind]
      cost.work <- cost.work.high
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuL-sH-wH-sdM'){
      cost.nu <- cost.nu.low
      cost.learning <- cost.data$cost.learning.high[cost.ind]
      cost.work <- cost.work.high
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuL-sL-wL-sdM'){
      cost.nu <- cost.nu.low
      cost.learning <- cost.data$cost.learning.low[cost.ind]
      cost.work <- cost.work.low
      cost.social <- cost.social.mid      
    }else if(path == 'oc-nuL-sH-wL-sdM'){
      cost.nu <- cost.nu.low
      cost.learning <- cost.data$cost.learning.high[cost.ind]
      cost.work <- cost.work.low
      cost.social <- cost.social.mid      
    }
    
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
      normal.school.year = tail(normal.school.year,n.weeks),
      n.school.weeks = sum(1-tail(normal.school.year,n.weeks))
    )    
    
    # load OC results
    u.noschool.obs <- list()
    u.noschool.full <- list()
    
    u.oc.tmp[[j]] <-
      as.matrix(
        read_csv(
          paste0('oc-results/',abbrev,'-',file.path,'-opt.csv'),
          show_col_types = FALSE)[,-1])    
    
    u.noschool.obs[[state.num]] <- params$u
    u.noschool.obs[[state.num]][1,] <- tail(normal.school.year,n.weeks)
    
    u.full <- matrix(1,nrow=df$p,ncol=n.weeks)
    u.noschool.full[[state.num]] <- u.full
    u.noschool.full[[state.num]][1,] <- tail(normal.school.year,n.weeks)      
    
    u.open <- matrix(0,nrow=df$p,ncol=n.weeks)
    u.open[1,] <- tail(normal.school.year,n.weeks)
    
    cost.oc[j,state.num] <- cost(u.oc.tmp[[j]],params)*params$n.pop
    cost.open[j,state.num] <- cost(u.open,params)*params$n.pop
    cost.obs[j,state.num] <- cost(params$u,params)*params$n.pop
    cost.full[j,state.num] <- cost(u.full,params)*params$n.pop
    cost.noschool.full[j,state.num] <- cost(u.noschool.full[[state.num]],params)*params$n.pop
    cost.noschool.obs[j,state.num] <- cost(u.noschool.obs[[state.num]],params)*params$n.pop
  }
  u.oc[[state.num]] <- u.oc.tmp
}

# plot costs
layout.matrix <- matrix(1:8, nrow = 4, ncol = 2,byrow=TRUE)
layout(mat = layout.matrix,
       heights = rep(1,4), # Heights of the rows
       widths = c(1, 1)) # Widths of the columns
par(mar=c(2,4,3,1))

titles <- c(
  'VSCD high, learning low, work high',
  'VSCD high, learning high, work high',
  'VSCD high, learning low, work low',
  'VSCD high, learning high, work low',
  'VSCD low, learning low, work high',
  'VSCD low, learning high, work high',
  'VSCD low, learning low, work low',
  'VSCD low, learning high, work low'  
)

for(j in 1:length(paths)){
  boxplot(
    cbind(
      log10(cost.oc[j,]/stan.data$n),
      log10(cost.full[j,]/stan.data$n),
      log10(cost.obs[j,]/stan.data$n),
      log10(cost.noschool.obs[j,]/stan.data$n),
      log10(cost.open[j,]/stan.data$n)
    ), 
    notch = TRUE,
    col=c(
      rgb(0,1,0,0.3),
      rgb(0,0,1,0.3),
      rgb(0,0,0,0.3),
      rgb(1,0,1,0.3),
      rgb(1,0,0,0.3)
    ), 
    ylab = TeX('$\\log_{10}(Cost)$ (\\$ per capita)'),
    main = titles[j],
    names = c('OC','Full','Obs.','Obs. - school','Open')
  )
}

# plot means
layout.matrix <- matrix(1:8, nrow = 4, ncol = 2,byrow=TRUE)
layout(mat = layout.matrix,
       heights = rep(1,4), # Heights of the rows
       widths = c(1, 1)) # Widths of the columns

titles <- c(
  'VSCD high, learning low, work high',
  'VSCD high, learning high, work high',
  'VSCD high, learning low, work low',
  'VSCD high, learning high, work low',
  'VSCD low, learning low, work high',
  'VSCD low, learning high, work high',
  'VSCD low, learning low, work low',
  'VSCD low, learning high, work low'  
)

u.oc.mean <- list()
for(j in 1:length(paths)){
  par(mar=c(2,4,3,1))
  
  u.oc.mean[[j]] <- 
    sapply(u.oc,function(x){rowMeans(x[[j]])}) 
  u.oc.mean[[j]][1,] <- sapply(
    u.oc,function(x){
      n.weeks <- length(x[[j]][1,])
      return(
        sum(x[[j]][1,] - tail(normal.school.year,n.weeks))/        
          (n.weeks - sum(tail(normal.school.year,n.weeks)))
      )
    }
  )
  ylab <- 'NPI average value'
  yaxt <- NULL
  if(j %% 2){
    ylab <- 'NPI average value'
    par(mar=c(2,4,3,1))
    yaxt <- NULL
  }
  
  boxplot(
    t(u.oc.mean[[j]])[,c(1,2,3,9,10,11)],  ylim = c(0,1),
    names = c('School','Work','Social','Test','Trace','Mask'),
    horizontal=FALSE,cex.axis=1, ylab = ylab, yaxt = yaxt,
    main = titles[j]
  )
}
