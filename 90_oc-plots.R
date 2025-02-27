rm(list=ls())

# plot optimal control and other policies

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

# generate results
u.obs <- list()
for(state in 1:51){
  u.obs[[state]] <- rbind(
    df$u[(df$t[state]+1):df$t[state+1],1], #school
    df$u[(df$t[state]+1):df$t[state+1],2], #work
    rowMeans(df$u[(df$t[state]+1):df$t[state+1],3:8]), #social
    df$u[(df$t[state]+1):df$t[state+1],9], #tests
    df$u[(df$t[state]+1):df$t[state+1],10], #trace
    df$u[(df$t[state]+1):df$t[state+1],11] #mask
    )
}

u.oc <- list()
u.open <- list()
u.noschool.obs <- list()
u.noschool.full <- list()

sim.oc <- list()
sim.obs <- list()
sim.open <- list()
sim.full <- list()
sim.noschool.obs <- list()
sim.noschool.full <- list()

cost.oc <- rep(NA,51)
cost.open <- rep(NA,51)
cost.obs <- rep(NA,51)
cost.full <- rep(NA,51)
cost.noschool.obs <- rep(NA,51)
cost.noschool.full <- rep(NA,51)

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

# cost of observed school closures in each state
cost.schools.obs <- rep(NA,51)

# states <- c(5,44,10,35) # CA, TX, FL, NY
states <- 1:51
nX <- 2
cost.setting <- 'nuL-sL-wM'
for(state.num in states){
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
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.mid
  cost.social <- cost.social.mid
  
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
  
  u.noschool.obs[[state.num]] <- params$u
  u.noschool.obs[[state.num]][1,] <- tail(normal.school.year,n.weeks)
  u.full <- matrix(1,nrow=df$p,ncol=n.weeks)
  u.noschool.full[[state.num]] <- u.full
  u.noschool.full[[state.num]][1,] <- tail(normal.school.year,n.weeks)  
  u.open[[state.num]] <- matrix(0,nrow=df$p,ncol=n.weeks)
  u.open[[state.num]][1,] <- tail(normal.school.year,n.weeks)
  
  sim.obs[[state.num]] <- list(
    d = seird.sim(params$u,params)$d,
    nu = seird.sim(params$u,params)$nu)
  sim.full[[state.num]] <- list(
    d = seird.sim(u.full,params)$d,
    nu = seird.sim(u.full,params)$nu)
  sim.open[[state.num]] <- list(
    d = seird.sim(u.open[[state.num]],params)$d,
    nu = seird.sim(u.open[[state.num]],params)$nu)
  sim.noschool.obs[[state.num]] <- list(
    d = seird.sim(u.noschool.obs[[state.num]],params)$d,
    nu = seird.sim(u.noschool.obs[[state.num]],params)$nu)
  sim.noschool.full[[state.num]] <- list(
    d = seird.sim(u.noschool.full[[state.num]],params)$d,
    nu = seird.sim(u.noschool.full[[state.num]],params)$nu)
  
  cost.open[state.num] <- cost.sim(u.open[[state.num]],sim.open[[state.num]]$nu,params)*params$n.pop
  cost.obs[state.num] <- cost.sim(params$u,sim.obs[[state.num]]$nu,params)*params$n.pop
  cost.full[state.num] <- cost.sim(u.full,sim.full[[state.num]]$nu,params)*params$n.pop
  cost.noschool.obs[state.num] <- cost.sim(u.noschool.obs[[state.num]],sim.noschool.obs[[state.num]]$nu,params)*params$n.pop
  cost.noschool.full[state.num] <- cost.sim(u.noschool.full[[state.num]],sim.noschool.full[[state.num]]$nu,params)*params$n.pop
  
  cost.schools.obs[state.num] <- cost.schools(sum(params$u[1,]),params)*params$n.pop
  
  u.oc[[state.num]] <-
    as.matrix(
      read_csv(
        paste0('oc-results/',abbrev,'-',path,'-opt.csv'),
        show_col_types = FALSE)[,-1])  
  
  sim.oc[[state.num]] <- list(
    d = seird.sim(u.oc[[state.num]],params)$d,
    nu = seird.sim(u.oc[[state.num]],params)$nu)
  cost.oc[state.num] <- cost.sim(u.oc[[state.num]],sim.oc[[state.num]]$nu,params)*params$n.pop
}

# plot OC results

par(mar=c(3,3,2,1.25))
layout.matrix <- matrix(c(1, 2, 3, 4,5,5), nrow = 3, ncol = 2,byrow=TRUE)

layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(1, 1)) # Widths of the two columns

# boxplot of deaths under each policy
boxplot(
  cbind(
    100*sapply(sim.oc,function(x){median(rowSums(x$d))})/stan.data$n,
    100*sapply(sim.full,function(x){median(rowSums(x$d))})/stan.data$n,
    100*sapply(sim.obs,function(x){median(rowSums(x$d))})/stan.data$n,
    100*sapply(sim.noschool.obs,function(x){median(rowSums(x$d))})/stan.data$n,
    100*sapply(sim.open,function(x){median(rowSums(x$d))})/stan.data$n
  ), ylim = c(0,0.8),
  notch = TRUE,
  col=c(
    rgb(0,1,0,0.3),
    rgb(0,0,1,0.3),
    rgb(0,0,0,0.3),
    rgb(1,0,1,0.3),
    rgb(1,0,0,0.3)
  ), main = 'Deaths per 100 pop.',
  names = c(
    'OC',
    'Full',
    'Obs.','Obs. - school','Open')
)

# boxplot of costs of each policy
boxplot(
  cbind(
    cost.oc/stan.data$n,
    cost.full/stan.data$n,
    cost.obs/stan.data$n,
    cost.noschool.obs/stan.data$n,
    cost.open/stan.data$n
  )/1000, ylim = c(0,40),
  notch = TRUE,
  col=c(
    rgb(0,1,0,0.3),
    rgb(0,0,1,0.3),
    rgb(0,0,0,0.3),
    rgb(1,0,1,0.3),
    rgb(1,0,0,0.3)
  ), main = 'Cost ($1000s per capita)',
  names = c(
    'OC',
    'Full',
    'Obs.','Obs. - school','Open')
)
u.oc.mean <- sapply(u.oc,rowMeans)
u.oc.mean[1,] <- sapply(
  u.oc,function(x){
    n.weeks <- length(x[1,])
    return(
      sum(x[1,] - tail(normal.school.year,n.weeks))/
        (n.weeks - sum(tail(normal.school.year,n.weeks)))
    )
  }
)

# boxplots of policy strength
boxplot(
  t(u.oc.mean)[,c(1,2,3,9,10,11)],  ylim = c(0,1),
  names = c('School','Work','Social','Test','Trace','Mask'),
  horizontal=FALSE,cex.axis=1,
  main = 'Average strength of NPI in optimal strategy'
)
u.obs.mean <- sapply(u.obs,rowMeans)
boxplot(
  t(u.obs.mean),  ylim = c(0,1),
  names = c('School','Work','Social','Test','Trace','Mask'),
  horizontal=FALSE,cex.axis=1,
  main = 'Average strength of NPI in observed strategy'
)

# plot of optimal workplace closure
u.oc.work <- sapply(u.oc,function(x){x[2,]})
n.weeks <- sapply(u.oc.work,length)
max.weeks <- 52
work <- matrix(0,nrow=51,ncol=max.weeks)
for(state in 1:51){
  work[state,(max.weeks-n.weeks[state]+1):max.weeks] <- u.oc.work[[state]]
}
work.int <- apply(work,2,function(x){quantile(x,probs=c(0.25,0.5,0.75))})
weeks <- lubridate::ymd( "2020-01-01" ) + lubridate::weeks(1:52)
plot(weeks,
     colMeans(work),
     type='l',main = 'Average strength of optimal workplace closure across states',
     xlab = 'Date',ylim=c(0,1))

#### plot trajectories of deaths under various policies

# aggregate deaths across states
n.days <- stan.data$n_d
max.days <- max(n.days)
dates <- ymd(model.data$dates[stan.data$state_inds==which.max(n.days)])
d.oc <- matrix(0,nrow=n.iter,ncol=max.days)
d.obs <- matrix(0,nrow=n.iter,ncol=max.days)
d.noschool.obs <- matrix(0,nrow=n.iter,ncol=max.days)
d.full <- matrix(0,nrow=n.iter,ncol=max.days)
d.noschool.full <- matrix(0,nrow=n.iter,ncol=max.days)
d.open <- matrix(0,nrow=n.iter,ncol=max.days)
d <- rep(0,max.days)
for(state in 1:51){
  d.oc[,(max.days-n.days[state]+1):max.days] <-
    d.oc[,(max.days-n.days[state]+1):max.days] +
    sim.oc[[state]]$d
  
  d.obs[,(max.days-n.days[state]+1):max.days] <- 
    d.obs[,(max.days-n.days[state]+1):max.days] +
    sim.obs[[state]]$d
  
  d.noschool.obs[,(max.days-n.days[state]+1):max.days] <- 
    d.noschool.obs[,(max.days-n.days[state]+1):max.days] +
    sim.noschool.obs[[state]]$d
  
  d.open[,(max.days-n.days[state]+1):max.days] <- 
    d.open[,(max.days-n.days[state]+1):max.days] +
    sim.open[[state]]$d
  
  d.full[,(max.days-n.days[state]+1):max.days] <- 
    d.full[,(max.days-n.days[state]+1):max.days] +
    sim.full[[state]]$d
  
  d.noschool.full[,(max.days-n.days[state]+1):max.days] <- 
    d.noschool.full[,(max.days-n.days[state]+1):max.days] +
    sim.noschool.full[[state]]$d
  
  d[(max.days-n.days[state]+1):max.days] <- 
    d[(max.days-n.days[state]+1):max.days] +
    stan.data$d[stan.data$state_inds == state]
}

quantile(rowSums(d.open),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.full),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.oc),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.obs),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.noschool.obs),probs=c(0.025,0.5,0.975))

quantile(rowSums(d.full)/rowSums(d.obs),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.noschool.obs)-rowSums(d.obs),probs=c(0.025,0.5,0.975))
quantile(2e6/(rowSums(d.noschool.obs)-rowSums(d.obs)),probs=c(0.025,0.5,0.975))
quantile(rowSums(d.obs)-rowSums(d.oc),probs=c(0.025,0.5,0.975))

# probability optimal policy saves lives
mean(rowSums(d.oc) < sum(d))
mean(rowSums(d.oc) < rowSums(d.obs))
mean(rowSums(d.oc) < sum(d) - 100000)
mean(rowSums(d.oc) < rowSums(d.obs) - 100000)

# ICER for Obs. - school vs. Obs.
(sum(cost.obs,na.rm=TRUE) - sum(cost.noschool.obs,na.rm=TRUE))/
  (1e6*(median(rowSums(d.noschool.obs)) - sum(d)))
(sum(cost.obs,na.rm=TRUE) - sum(cost.noschool.obs,na.rm=TRUE))/
  (1e6*(mean(rowSums(d.noschool.obs)) - sum(d)))
(sum(cost.obs,na.rm=TRUE) - sum(cost.noschool.obs,na.rm=TRUE))/
  (1e6*(median(rowSums(d.noschool.obs) - rowSums(d.obs))))
(sum(cost.obs,na.rm=TRUE) - sum(cost.noschool.obs,na.rm=TRUE))/
  (1e6*(mean(rowSums(d.noschool.obs) - rowSums(d.obs))))
quantile((sum(cost.obs,na.rm=TRUE) - sum(cost.noschool.obs,na.rm=TRUE))/
  (1e6*(rowSums(d.noschool.obs) - rowSums(d.obs))),
  probs = c(0.5,0.025,0.975))

# make the plot
par(mfrow=c(2,1))
par(mar=c(0,4,5,1))
p <- 0.05
d.obs.hat <- matrix(
  apply(
    d.obs,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.noschool.obs.hat <- matrix(
  apply(
    d.noschool.obs,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.oc.hat <- matrix(
  apply(
    d.oc,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.full.hat <- matrix(
  apply(
    d.full,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.noschool.full.hat <- matrix(
  apply(
    d.noschool.full,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.open.hat <- matrix(
  apply(
    d.open,
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
plot(dates,
     log10(d.obs.hat)[2,], 
     type='l', 
     ylim = c(1,0.5+max(
       log10(d.obs.hat),
       log10(d.noschool.obs.hat),
       log10(d.full.hat),
       log10(d.noschool.full.hat),
       log10(d.oc.hat),
       log10(d.open.hat))),
     xlab = 'Month (2020)',
     xaxt = 'n',
     ylab = TeX('$\\log_{10}$(Daily deaths)'),
     xlim=as.Date(c("2020-03-08", "2020-12-31"))
)
mtext('Deaths under various policies',font=2, side = 3, line = -4, outer = TRUE)
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.obs.hat)[1,]),log10(d.obs.hat)[3,]),
  border = NA,
  col = rgb(0,0,0,0.2)
)
points(dates,log10(d),pch='.')
lines(dates,
      log10(d.noschool.obs.hat)[2,],
      col=rgb(1,0,1,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.noschool.obs.hat)[1,]),log10(d.noschool.obs.hat)[3,]),
  border = NA,
  col = rgb(1,0,1,0.2)
)
lines(dates,
      log10(d.open.hat)[2,],
      col=rgb(1,0,0,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.open.hat)[1,]),log10(d.open.hat)[3,]),
  border = NA,
  col = rgb(1,0,0,0.2)
)
lines(dates,
      log10(d.full.hat)[2,],
      col=rgb(0,0,1,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.full.hat)[1,]),log10(d.full.hat)[3,]),
  border = NA,
  col = rgb(0,0,1,0.2)
)
lines(dates,
      log10(d.oc.hat)[2,],
      col=rgb(0,1,0,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.oc.hat)[1,]),log10(d.oc.hat)[3,]),
  border = NA,
  col = rgb(0,1,0,0.2)
)

d.obs.hat <- matrix(
  apply(
    t(apply(d.obs,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.noschool.obs.hat <- matrix(
  apply(
    t(apply(d.noschool.obs,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.oc.hat <- matrix(
  apply(
    t(apply(d.oc,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.full.hat <- matrix(
  apply(
    t(apply(d.full,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.noschool.full.hat <- matrix(
  apply(
    t(apply(d.noschool.full,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
d.open.hat <- matrix(
  apply(
    t(apply(d.open,1,cumsum)),
    2,function(x){quantile(x,probs=c(p,0.5,1-p))}),
  nrow = 3, byrow=FALSE
)
par(mar=c(4,4,1,1))
plot(dates,
     log10(d.obs.hat)[2,], 
     type='l', 
     ylim = c(1,0.5+max(
       log10(d.obs.hat),
       log10(d.noschool.obs.hat),
       log10(d.full.hat),
       log10(d.noschool.full.hat),
       log10(d.oc.hat),
       log10(d.open.hat))),
     xlab = 'Month (2020)', 
     ylab = TeX('$\\log_{10}$(Cumulative deaths)'),
     xlim=as.Date(c("2020-03-08", "2020-12-31"))
)
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.obs.hat)[1,]),log10(d.obs.hat)[3,]),
  border = NA,
  col = rgb(0,0,0,0.2)
)
points(dates,log10(cumsum(d)),pch='.') 
lines(dates,
      log10(d.noschool.obs.hat)[2,],
      col=rgb(1,0,1,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.noschool.obs.hat)[1,]),log10(d.noschool.obs.hat)[3,]),
  border = NA,
  col = rgb(1,0,1,0.2)
)
lines(dates,
      log10(d.open.hat)[2,],
      col=rgb(1,0,0,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.open.hat)[1,]),log10(d.open.hat)[3,]),
  border = NA,
  col = rgb(1,0,0,0.2)
)
lines(dates,
      log10(d.full.hat)[2,],
      col=rgb(0,0,1,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.full.hat)[1,]),log10(d.full.hat)[3,]),
  border = NA,
  col = rgb(0,0,1,0.2)
)
lines(dates,
      log10(d.oc.hat)[2,],
      col=rgb(0,1,0,1))
polygon(
  x = c(rev(dates),dates),
  y = c(rev(log10(d.oc.hat)[1,]),log10(d.oc.hat)[3,]),
  border = NA,
  col = rgb(0,1,0,0.2)
)
legend(
  'bottomright',
  legend=c(
    paste0('Open, D = ',
           round(median(rowSums(d.open))),
           round(sum(cost.open,na.rm=TRUE)/sum(stan.data$n))),
    paste0('Obs. - school, D = ',
           round(median(rowSums(d.noschool.obs))),
           round(sum(cost.noschool.obs,na.rm=TRUE)/sum(stan.data$n))),
    paste0('Observed, D = ',
           round(median(rowSums(d.obs))),           
           round(sum(cost.obs,na.rm=TRUE)/sum(stan.data$n))),
    paste0('Optimal, D = ',
           round(median(rowSums(d.oc))),
           round(sum(cost.oc)/sum(stan.data$n))),
    paste0('Full lockdown, D = ',
           round(median(rowSums(d.full))),
           round(sum(cost.full,na.rm=TRUE)/sum(stan.data$n)))
  ), 
  col = c(
    'red',
    'purple',
    'black',    
    'green',
    'blue'
  ),cex = 0.85,lty=1)




##### calculate ICERs over time
states <- c(5,44,10,35) # CA, TX, FL, NY
icer.obs <- list()
icer.noschool.obs <- list()
icer.open <- list()
icer.full <- list()
icer.oc <- list()

for(state.num in states){
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
  cost.nu <- cost.nu.low
  cost.learning <- cost.data$cost.learning.low[cost.ind]
  cost.work <- cost.work.mid
  cost.social <- cost.social.mid
  
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
  
  cutoff <- 1 # start plotting after (cutoff) weeks
  icer.obs.tmp <- list()
  icer.noschool.obs.tmp <- list()
  icer.open.tmp <- list()
  icer.full.tmp <- list()
  icer.oc.tmp <- list()
  
  u.full <- matrix(1,nrow=df$p,ncol=n.weeks)
  for(j in (1+cutoff):n.weeks){
    print(100*j/n.weeks)
    
    # calculate infections averted (compared to baseline fully open policy)
    icer.open.tmp[[j-cutoff]] <- cost.cumulative.sim(
      u.open[[state.num]],sim.open[[state.num]]$nu,params,j)
    
    # optimal control
    icer.oc.tmp[[j-cutoff]] <- cost.cumulative.sim(
      u.oc[[state.num]],sim.oc[[state.num]]$nu,params,j)
    
    # observed
    icer.obs.tmp[[j-cutoff]] <- cost.cumulative.sim(
      params$u,sim.obs[[state.num]]$nu,params,j)
    
    # observed - schools
    icer.noschool.obs.tmp[[j-cutoff]] <- cost.cumulative.sim(
      u.noschool.obs[[state.num]],sim.noschool.obs[[state.num]]$nu,params,j)
    
    # full lockdown
    icer.full.tmp[[j-cutoff]] <- cost.cumulative.sim(
      u.full,sim.full[[state.num]]$nu,params,j)
  }
  
  icer.oc[[state.num]] <- icer.oc.tmp
  icer.full[[state.num]] <- icer.full.tmp
  icer.open[[state.num]] <- icer.open.tmp
  icer.obs[[state.num]] <- icer.obs.tmp
  icer.noschool.obs[[state.num]] <- icer.noschool.obs.tmp
}

cumulative <- FALSE

par(mfrow=c(2,2))
par(mar=c(4,4,4,1))

icer.baseline <- icer.open
for(state in states){
  n.weeks <- length(icer.obs[[state]])
  
  icer.obs.tmp <- list()
  icer.noschool.obs.tmp <- list()
  icer.oc.tmp <- list()
  icer.full.tmp <- list()
  icer.open.tmp <- list()
  
  for(j in 1:n.weeks){
    if(j==1 || cumulative){
      icer.obs.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.obs[[state]][[j]]$cost.inf) /
        (icer.obs[[state]][[j]]$cost.npis - icer.baseline[[state]][[j]]$cost.npis)
      
      icer.noschool.obs.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.noschool.obs[[state]][[j]]$cost.inf) /
        (icer.noschool.obs[[state]][[j]]$cost.npis - icer.baseline[[state]][[j]]$cost.npis)
      
      icer.oc.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.oc[[state]][[j]]$cost.inf) /
        (icer.oc[[state]][[j]]$cost.npis - icer.baseline[[state]][[j]]$cost.npis)
      
      icer.open.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.open[[state]][[j]]$cost.inf) /
        (icer.open[[state]][[j]]$cost.npis - icer.baseline[[state]][[j]]$cost.npis)
      
      icer.full.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.full[[state]][[j]]$cost.inf) /
        (icer.full[[state]][[j]]$cost.npis - icer.baseline[[state]][[j]]$cost.npis)
    }else{
      icer.obs.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.baseline[[state]][[j-1]]$cost.inf  
         - (icer.obs[[state]][[j]]$cost.inf - icer.obs[[state]][[j-1]]$cost.inf)) /
        (icer.obs[[state]][[j]]$cost.npis - icer.obs[[state]][[j-1]]$cost.npis  
         - (icer.baseline[[state]][[j]]$cost.npis - icer.baseline[[state]][[j-1]]$cost.npis))
      
      icer.noschool.obs.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.baseline[[state]][[j-1]]$cost.inf  
         - (icer.noschool.obs[[state]][[j]]$cost.inf - icer.noschool.obs[[state]][[j-1]]$cost.inf)) /
        (icer.noschool.obs[[state]][[j]]$cost.npis - icer.noschool.obs[[state]][[j-1]]$cost.npis  
         - (icer.baseline[[state]][[j]]$cost.npis - icer.baseline[[state]][[j-1]]$cost.npis))  
      
      icer.full.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.baseline[[state]][[j-1]]$cost.inf  
         - (icer.full[[state]][[j]]$cost.inf - icer.full[[state]][[j-1]]$cost.inf)) /
        (icer.full[[state]][[j]]$cost.npis - icer.full[[state]][[j-1]]$cost.npis  
         - (icer.baseline[[state]][[j]]$cost.npis - icer.baseline[[state]][[j-1]]$cost.npis)) 
      
      icer.oc.tmp[[j]] <- 
        (icer.baseline[[state]][[j]]$cost.inf - icer.baseline[[state]][[j-1]]$cost.inf  
         - (icer.oc[[state]][[j]]$cost.inf - icer.oc[[state]][[j-1]]$cost.inf)) /
        (icer.oc[[state]][[j]]$cost.npis - icer.oc[[state]][[j-1]]$cost.npis  
         - (icer.baseline[[state]][[j]]$cost.npis - icer.baseline[[state]][[j-1]]$cost.npis))
      
      icer.open.tmp[[j]] <-
        (icer.baseline[[state]][[j]]$cost.inf - icer.baseline[[state]][[j-1]]$cost.inf
         - (icer.open[[state]][[j]]$cost.inf - icer.open[[state]][[j-1]]$cost.inf)) /
        (icer.open[[state]][[j]]$cost.npis - icer.open[[state]][[j-1]]$cost.npis
         - (icer.baseline[[state]][[j]]$cost.npis - icer.baseline[[state]][[j-1]]$cost.npis))
    }    
  }
  
  weeks <- rev(lubridate::ymd( "2020-12-31" ) - lubridate::weeks(1:n.weeks))
  
  p <- 0.25
  icer.obs.hat <- 
    matrix(unlist(lapply(icer.obs.tmp,function(x){
      quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)})),
      nrow = 3, byrow=FALSE)
  icer.oc.hat <- 
    matrix(
      unlist(lapply(icer.oc.tmp,function(x){
        quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)})),
      nrow = 3, byrow=FALSE
    )
  icer.noschool.obs.hat <- 
    matrix(
      unlist(lapply(icer.noschool.obs.tmp,function(x){
        quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)})),
      nrow = 3, byrow=FALSE
    )
  icer.full.hat <- 
    matrix(
      unlist(lapply(icer.full.tmp,function(x){
        quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)})),
      nrow = 3, byrow=FALSE
    )
  icer.open.hat <- 
    matrix(
      unlist(lapply(icer.open.tmp,function(x){
        quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)})),
      nrow = 3, byrow=FALSE
    )
  
  if((state == 44) & cumulative){#TX
    plot(
      weeks,icer.obs.hat[2,],type='l', 
      ylim = c(0,15),
      xlab = 'Month (2020)',
      xlim = ymd(c('2020-03-19','2020-12-24')),
      ylab = 'Standardized ICER',
      main = model.data$states[state]
    )      
  }else if(state == 44){#TX  
    plot(
      weeks,icer.obs.hat[2,],type='l', 
      ylim = c(0,40),
      xlab = 'Month (2020)',
      xlim = ymd(c('2020-03-19','2020-12-24')),
      ylab = 'Standardized ICER',
      main = model.data$states[state]
    )    
  }else{
    plot(
      weeks,icer.obs.hat[2,],type='l', 
      ylim = c(
        min(0,
            icer.oc.hat,
            icer.obs.hat,
            icer.noschool.obs.hat,
            icer.full.hat)
        ,
        max(
          icer.oc.hat,
          icer.obs.hat,
          icer.noschool.obs.hat,
          icer.full.hat,
          na.rm=TRUE)),
      xlab = 'Month (2020)',
      ylab = 'Standardized ICER',
      main = model.data$states[state]
    )      
    }
  abline(h = 1, lty=2,col='grey')
  polygon(
    x = c(rev(weeks),weeks),
    y = c(rev(icer.obs.hat[1,]),icer.obs.hat[3,]),
    border = NA,
    col = rgb(0,0,0,0.2)
  )
  lines(weeks,icer.noschool.obs.hat[2,],col='purple')
  polygon(
    x = c(rev(weeks),weeks),
    y = c(rev(icer.noschool.obs.hat[1,]),icer.noschool.obs.hat[3,]),
    border = NA,
    col = rgb(1,0,1,0.2)
  )
  lines(weeks,icer.full.hat[2,],col='blue')
  polygon(
    x = c(rev(weeks),weeks),
    y = c(rev(icer.full.hat[1,]),icer.full.hat[3,]),
    border = NA,
    col = rgb(0,0,1,0.2)
  )
  lines(weeks,icer.oc.hat[2,],col='green')
  polygon(
    x = c(rev(weeks),weeks),
    y = c(rev(icer.oc.hat[1,]),icer.oc.hat[3,]),
    border = NA,
    col = rgb(0,1,0,0.2)
  )
}

legend(
  'topright',
  legend = c('OC','Obs. - school','Obs.','Full'),
  lty = 1,
  col = c('green','purple','grey','blue')
)

