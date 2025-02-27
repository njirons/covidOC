rm(list=ls())

# posterior plots

library(dplyr)
library(rjson)
library(cmdstanr)
library(lubridate)
library(vioplot)
library(latex2exp)

# read in data
stan_data <- readRDS(file="data/stan_data.RData")
stan.data <- stan_data
model_data <- readRDS(file="data/model_data.RData")
model.data <- model_data

# violinplot of NPI effects
par(mfrow=c(1,1),
    mar=c(2,4,0.5,0.5))
npi.names <- c(
  'Schools',
  'Work',
  'Events',
  'Gatherings',
  'Transit',
  'Stay-at-home',
  'Movement',
  'Public Info.',
  'Testing',
  'Tracing',
  'Masking'
)
load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',2,'-full.Rda'))
mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
inds <- order(apply(mu,2,median))
vioplot(
  100*(1-exp(mu[,inds])),
  names = npi.names[inds],
  cex.names = 1,
  ylim=c(0,40),
  ylab = TeX('Percent reduction in $R_0$')
  )
abline(h=10,lty=2,
       col=rgb(0,0,0,0.3))
abline(h=20,lty=2,
       col=rgb(0,0,0,0.3)) 
abline(h=30,lty=2,
       col=rgb(0,0,0,0.3)) 

par(mfrow=c(3,1))
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  inds <- order(apply(mu,2,median))
  if(j==1){
    par(mar=c(2,2,0.5,0.5))
    vioplot(
      100*(1-exp(mu[,inds])),
      names = npi.names[inds],
      cex.names = 1.0,
      ylim=c(0,40),
      ylab = TeX('Percent reduction in $R_0$')
      ) 
    abline(h=10,lty=2,
           col=rgb(0,0,0,0.3))
    abline(h=20,lty=2,
           col=rgb(0,0,0,0.3)) 
    abline(h=30,lty=2,
           col=rgb(0,0,0,0.3))     
  }
  else{
    par(mar=c(2,2,0.5,0.5))
    vioplot(
      100*(1-exp(mu[,inds])),
      names = npi.names[inds],
      cex.names = 1.0,
      ylim=c(0,40),
      ylab = TeX('Percent reduction in $R_0$')
      )
    abline(h=10,lty=2,
           col=rgb(0,0,0,0.3))
    abline(h=20,lty=2,
           col=rgb(0,0,0,0.3)) 
    abline(h=30,lty=2,
           col=rgb(0,0,0,0.3))    
  }
}

# effects of individual NPIs
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  inds <- order(apply(mu,2,median))  
  print(apply(100*(1-exp(mu[,inds])),2,function(x){quantile(x,c(0.5,0.025,0.975))}))
}

# effects of individual NPIs
ind <- 10
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- as.matrix((sims.npi %>% select(starts_with('mu[')))[,2:12])
  print(quantile(100*(1-exp(mu[,ind])),probs=c(0.5,0.025,0.975)))
}

# effects of combined social distancing measures
inds <- 3:8
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  print(quantile(100*(1-exp(rowSums(mu[,inds]))),probs=c(0.5,0.025,0.975)))
}

# effects of deaths/removals/infections
mu.dri <- (sims.npi %>% select(starts_with('mu[')))[,13:14]
vioplot(
  100*(1-exp(mu.dri)),
  # names = c('Deaths','Removals','Infections'),
  names = c('Deaths','Removals'),
  cex.names = 0.7,
  # ylim=c(0,40),
  main = TeX('Percent reduction in $R_0$ by endogenous response'))


# calculate national ifr
n.iter <- 4000
ifr.state <- list()
d.state <- rep(NA,51)
d.state.post1 <- matrix(NA,nrow=n.iter,ncol=51)
d.state.post2 <- matrix(NA,nrow=n.iter,ncol=51)
for(state in 1:51){
  print(state)
  abbrev <- model_data$abbrev[state]
  load(paste0('state-sims/sims_',abbrev,'-',state,'.Rda'))
  ifr.state[[state]] <- sims$ifr
  d.state[state] <- stan_data$d0[state]+
    sum(stan_data$d[(stan_data$t_d[state]+1):stan_data$t_d[state+1]])
  sus <- as.matrix(sims %>% select(starts_with('sus[')))
  exps <- as.matrix(sims %>% select(starts_with('exps[')))
  inf <- as.matrix(sims %>% select(starts_with('inf[')))
  rem_s <- as.matrix(sims %>% select(starts_with('rem_s[')))
  rem_d <- as.matrix(sims %>% select(starts_with('rem_d[')))
  ind <- dim(sus)[2]
  
  d.state.post1[,state] <- stan_data$n[state] * (
    1 - (
      sus[,ind] + exps[,ind] + inf[,ind] + rem_s[,ind]
    )
  )
  
  d.state.post2[,state] <- stan_data$n[state] * (
    1 - (
      sus[,ind] + exps[,ind] + inf[,ind] + rem_s[,ind] + rem_d[,ind]
    )
  )  
}
sum(d.state)

ifr.usa <- 0
ifr.usa.post1 <- 0
ifr.usa.post2 <- 0
for(state in 1:51){
  ifr.usa <- ifr.usa + ifr.state[[state]]*d.state[state]
  ifr.usa.post1 <- ifr.usa.post1 + ifr.state[[state]]*d.state.post1[,state]
  ifr.usa.post2 <- ifr.usa.post2 + ifr.state[[state]]*d.state.post2[,state]
}
ifr.usa <- ifr.usa/sum(d.state)
ifr.usa.post1 <- ifr.usa.post1/rowSums(d.state.post1)
ifr.usa.post2 <- ifr.usa.post2/rowSums(d.state.post2)

round(quantile(100*ifr.usa,c(0.5,0.025,0.975)),2)
round(quantile(100*ifr.usa.post1,c(0.5,0.025,0.975)),2)
round(quantile(100*ifr.usa.post2,c(0.5,0.025,0.975)),2)

# calculate national infection/infectious/cases/deaths counts
n.iter <- 4000
n.days <- stan_data$n_d
max.days <- max(n.days)
dates <- ymd(model_data$dates[stan_data$state_inds==which.max(n.days)])
nu <- matrix(0,ncol=max.days,nrow=n.iter)
inf <- matrix(0,ncol=max.days,nrow=n.iter)
d <- rep(0,max.days)
d.par <- matrix(0,ncol=max.days,nrow=n.iter)
d.hat <- matrix(0,ncol=max.days,nrow=n.iter)
c <- rep(0,max.days)
c.par <- matrix(0,ncol=max.days,nrow=n.iter)
c.hat <- matrix(0,ncol=max.days,nrow=n.iter)
sus <- matrix(0,ncol=max.days,nrow=n.iter)
for(state in 1:51){
  print(state)
  abbrev <- model_data$abbrev[state]
  load(paste0('state-sims/sims_',abbrev,'-',state,'.Rda'))  
  
  d[(max.days-n.days[state]+1):max.days] <- 
    d[(max.days-n.days[state]+1):max.days] +
    stan_data$d[stan_data$state_inds == state]

  d.par[,(max.days-n.days[state]+1):max.days] <- 
    d.par[,(max.days-n.days[state]+1):max.days] + 
    as.matrix(sims %>% select(starts_with('death_par[')))
  
  d.hat[,(max.days-n.days[state]+1):max.days] <- 
    d.hat[,(max.days-n.days[state]+1):max.days] + 
    as.matrix(sims %>% select(starts_with('dhat_nz['))) * 
    matrix(
      sapply(sims$d_hurdle,function(x){rbinom(n.days[state],1,1-x)}),
      ncol = n.days[state])
  
  c[(max.days-n.days[state]+1):max.days] <- 
    c[(max.days-n.days[state]+1):max.days] +
    stan_data$c[stan_data$state_inds == state]
  
  c.par[,(max.days-n.days[state]+1):max.days] <- 
    c.par[,(max.days-n.days[state]+1):max.days] + 
    as.matrix(sims %>% select(starts_with('case_par[')))
  
  c.hat[,(max.days-n.days[state]+1):max.days] <- 
    c.hat[,(max.days-n.days[state]+1):max.days] + 
    as.matrix(sims %>% select(starts_with('chat_nz['))) * 
    matrix(
      sapply(sims$c_hurdle,function(x){rbinom(n.days[state],1,1-x)}),
      ncol = n.days[state])  
  
  inf[,(max.days-n.days[state]+1):max.days] <- 
    inf[,(max.days-n.days[state]+1):max.days] + 
    as.matrix((sims %>% select(starts_with('inf[')))*stan_data$n[state])  
  
  nu[,(max.days-n.days[state]+1):max.days] <- 
    nu[,(max.days-n.days[state]+1):max.days] + 
    as.matrix((sims %>% select(starts_with('nu[')))*stan_data$n[state])  
  
  sus[,(max.days-n.days[state]+1):max.days] <- 
    sus[,(max.days-n.days[state]+1):max.days] + 
    as.matrix((sims %>% select(starts_with('sus[')))*stan_data$n[state])
}

# total number of infections
quantile(rowSums(nu),probs=c(0.5,0.025,0.975))
100*quantile(rowSums(nu),probs=c(0.5,0.025,0.975))/(329.5*1e6)
100*quantile(rowSums(nu),probs=c(0.5,0.025,0.975))/sum(stan_data$n)
100*(1-quantile(sus[,max.days],probs=c(0.5,0.975,0.025))/(329.5*1e6))
100*(1-quantile(sus[,max.days],probs=c(0.5,0.975,0.025))/sum(stan_data$n))
100*(sum(stan_data$n)-quantile(sus[,max.days],probs=c(0.5,0.975,0.025)))

# calculate posterior R^2 of NPI regression model fit
bhm_data <- readRDS("data/bhm_data.RData")
empirical.variance <- 
  rep(unlist(lapply(bhm_data,function(x){var(x$log_rt)})),each=100)
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  eps <- sims.npi %>% select(contains('eps'))
  residual.variance <- apply(eps,1,var)
  print(summary(1-residual.variance/empirical.variance))  
}

# calculate posterior of AR(1) parameter from regression model
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  print(quantile(sims.npi$phi_ar,probs=c(0.5,0.025,0.975)))
}

# calculate posterior of degrees of freedom
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  print(quantile(exp(sims.npi$log_nu),probs=c(0.5,0.025,0.975)))
}

# calculate posterior of R0
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  print(quantile(sims.npi$R0,probs=c(0.5,0.025,0.975)))
}

# calculate total effect of NPIs
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  print(quantile(100*(1-exp(rowSums(mu))),probs=c(0.5,0.025,0.975)))
  print(quantile(sims.npi$R0*exp(rowSums(mu)),probs=c(0.5,0.025,0.975)))
}

# effects of deaths, recoveries, infections
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  print(quantile(sims.npi$`mu[13]`,probs=c(0.5,0.025,0.975)))
  if(j == 2){
    print(quantile(sims.npi$`mu[14]`,probs=c(0.5,0.025,0.975)))
  }
  if(j == 3){
    print(quantile(sims.npi$`mu[14]`,probs=c(0.5,0.025,0.975)))
    print(quantile(sims.npi$`mu[15]`,probs=c(0.5,0.025,0.975)))
  }
}

# calculate mean weekly death/removal/infection rate across states
d.mean <- mean(unlist(lapply(bhm_data,function(x){mean(x$X_all[,1])})))
sd(unlist(lapply(bhm_data,function(x){mean(x$X_all[,1])})))
r.mean <- mean(unlist(lapply(bhm_data,function(x){mean(x$X_all[,2])})))
sd(unlist(lapply(bhm_data,function(x){mean(x$X_all[,2])})))
i.mean <- mean(unlist(lapply(bhm_data,function(x){mean(x$X_all[,3])})))
sd(unlist(lapply(bhm_data,function(x){mean(x$X_all[,3])})))

# mean percent effects of deaths, recoveries, infections
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  if(j == 1){
    print(quantile(100*(1-exp(sims.npi$`mu[13]`*d.mean)),probs=c(0.5,0.025,0.975)))    
  }
  if(j == 2){
    print(
      quantile(100*(1-exp(
        sims.npi$`mu[13]`*d.mean +
          sims.npi$`mu[14]`*r.mean
        )),probs=c(0.5,0.025,0.975)))
  }
  if(j == 3){
    print(quantile(100*(1-exp(
      sims.npi$`mu[13]`*d.mean +
        sims.npi$`mu[14]`*r.mean +
        sims.npi$`mu[15]`*i.mean
      )),probs=c(0.5,0.025,0.975)))
  }
}

# effect of endogenous response as fraction of NPI effect
for(j in 1:3){
  load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  if(j == 1){
    print(quantile(
      (1-exp(sims.npi$`mu[13]`*d.mean))/
        (1-exp(rowSums(mu)))
      ,probs=c(0.5,0.025,0.975)))    
  }
  if(j == 2){
    print(
      quantile(
        (1-exp(sims.npi$`mu[13]`*d.mean +
          sims.npi$`mu[14]`*r.mean))/
          (1-exp(rowSums(mu)))
          ,probs=c(0.5,0.025,0.975)))
  }
  if(j == 3){
    print(quantile(
      (1-exp(sims.npi$`mu[13]`*d.mean +
        sims.npi$`mu[14]`*r.mean +
        sims.npi$`mu[15]`*i.mean))/
        (1-exp(rowSums(mu)))
      ,probs=c(0.5,0.025,0.975)))
  }
} 

# combined effect of endogenous response and NPIs
for(j in 1:3){
  load(paste0('npi-bhm-ar1-nX',j,'-full.Rda'))
  mu <- (sims.npi %>% select(starts_with('mu[')))[,2:12]
  if(j == 1){
    print(quantile(
      100*(1-exp(rowSums(mu) + sims.npi$`mu[13]`*d.mean))
      ,probs=c(0.5,0.025,0.975)))    
  }
  if(j == 2){
    print(quantile(
      100*(1-exp(
        rowSums(mu) + 
          sims.npi$`mu[13]`*d.mean + 
          sims.npi$`mu[14]`*r.mean))
      ,probs=c(0.5,0.025,0.975)))
  }
  if(j == 3){
    print(quantile(
      100*(1-exp(
        rowSums(mu) + 
          sims.npi$`mu[13]`*d.mean + 
          sims.npi$`mu[14]`*r.mean +
          sims.npi$`mu[15]`*i.mean))
      ,probs=c(0.5,0.025,0.975)))
  }
}     


# plot national deaths/infections/cases estimates
d.hat.ci <- apply(
  d.hat, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

d.par.ci <- apply(
  d.par, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

c.hat.ci <- apply(
  c.hat, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

c.par.ci <- apply(
  c.par, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

inf.ci <- apply(
  inf, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

nu.ci <- apply(
  nu, 2, function(x){quantile(x,probs=c(0.5,0.05,0.95))}
)

par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(
  dates,d/1e3,pch='.',
  xlab = 'Month', ylab = 'Deaths (thousands) / Infections (millions)'
)
polygon(c(rev(dates), dates), 
        c(rev(d.hat.ci[2,]), d.hat.ci[3,])/1e3,
        col = rgb(0,0,1,0.3), border = NA)
polygon(c(rev(dates), dates), 
        c(rev(d.par.ci[2,]), d.par.ci[3,])/1e3,
        col = rgb(0,0,1,0.5), border = NA)
lines(dates,d.par.ci[1,]/1e3,col='black')
polygon(c(rev(dates), dates), 
        c(rev(inf.ci[2,]), inf.ci[3,])/1e6,
        col = rgb(1,0,0,0.3), border = NA)
lines(dates,inf.ci[1,]/1e6,col='red')

plot(
  dates,c/1e3,pch='.',
  xlab = 'Month', ylab = 'Cases (thousands)'
)
polygon(c(rev(dates), dates), 
        c(rev(c.hat.ci[2,]), c.hat.ci[3,])/1e3,
        col = rgb(0,0,1,0.3), border = NA)
polygon(c(rev(dates), dates), 
        c(rev(c.par.ci[2,]), c.par.ci[3,])/1e3,
        col = rgb(0,0,1,0.5), border = NA)
lines(dates,c.par.ci[1,]/1e3,col='black')


### state-specific estimates
state_num <- 1

dates <- ymd(model_data$dates[stan_data$state_inds==state_num])
abbrev <- model_data$abbrev[state_num]
load(paste0('state-data/data_',abbrev,'-',state_num,'.Rda'))
load(paste0('state-sims/sims_',abbrev,'-',state_num,'.Rda'))

par(mfrow=c(2,2))
par(mar=c(0,4,5,0))

p <- 0.05
d_par <- apply(
  sims %>% select(starts_with('death_par[')),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

d_par_obs <- apply(
  sims %>% select(starts_with('death_par_obs[')),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

dhat <- apply(
  sims %>% select(starts_with('dhat_nz['))*
    rbinom(length(sims$d_hurdle),1,1-sims$d_hurdle),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

dhat_obs <- apply(
  (sims %>% select(starts_with('dhat_nz_obs[')))*
    rbinom(length(sims$d_hurdle),1,1-sims$d_hurdle),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

plot(dates,dhat[2,],ylim=range(state_data$d),type='l',col='red',ylab='Deaths',
     xaxt='n') 
polygon(c(rev(dates), dates), 
        c(rev(dhat_obs[1,]), dhat_obs[3,]),
        col = rgb(1,0,0,0.3), border = NA)
lines(dates,d_par[2,],col='blue')
polygon(c(rev(dates), dates),
        c(rev(d_par[1,]), d_par[3,]),
        col = rgb(0,0,1,0.3), border = NA)
points(dates,state_data$d,pch='.')

c_par <- apply(
  sims %>% select(starts_with('case_par[')),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

c_par_obs <- apply(
  sims %>% select(starts_with('case_par_obs[')),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

chat <- apply(
  sims %>% select(starts_with('chat_nz['))*
    rbinom(length(sims$c_hurdle),1,1-sims$c_hurdle),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

chat_obs <- apply(
  sims %>% select(starts_with('chat_nz_obs['))*
    rbinom(length(sims$c_hurdle),1,1-sims$c_hurdle),2,
  function(x){quantile(x,probs=c(p,0.5,1-p),na.rm=TRUE)}
)

plot(dates,c_par[2,],ylim=range(state_data$c),type='l',col='red',ylab='Cases',
     xaxt = 'n') 
polygon(c(rev(dates), dates), 
        c(rev(chat_obs[1,]), chat_obs[3,]), 
        col = rgb(1,0,0,0.3), border = NA)
polygon(c(rev(dates), dates),
        c(rev(c_par[1,]), c_par[3,]),
        col = rgb(0,0,1,0.3), border = NA)
points(dates,state_data$c,pch='.')

rt <- apply(
  sims %>% select(starts_with('r0[')),2,
  function(x){quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95))}
)

par(mar=c(4,4,1,0))
plot(dates,rt[3,state_data$error_term_inds],ylim=range(rt),type='l',col='black',
     ylab=TeX('$R_0(t)$'),xlab = 'Month (2020)')  
abline(h=1,col='grey',lty='dashed')
polygon(c(rev(dates), dates), 
        c(rev(rt[1,state_data$error_term_inds]), 
          rt[5,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.3), border = NA)
polygon(c(rev(dates), dates), 
        c(rev(rt[2,state_data$error_term_inds]), 
          rt[4,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.5), border = NA)
lines(dates,rt[3,state_data$error_term_inds])  

car <- apply(
  sims %>% select(starts_with('car[')),2,
  function(x){quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95))}
)

plot(dates,car[3,state_data$error_term_inds],
     ylim=c(0,1),type='l',col='black',xlab='Month (2020)',
     ylab = 'Case ascertainment rate')    
polygon(c(rev(dates), dates), 
        c(rev(car[1,state_data$error_term_inds]), 
          car[5,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.3), border = NA)
polygon(c(rev(dates), dates), 
        c(rev(car[2,state_data$error_term_inds]), 
          car[4,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.5), border = NA)
lines(dates,car[3,state_data$error_term_inds])  

mtext(model_data$states[state_num], 
      side = 3, line = -4, outer = TRUE,font=2)


inf <- apply(
  sims %>% select(starts_with('inf[')),2,
  function(x){quantile(x,probs=c(p,0.5,1-p))}
)
plot(dates,inf[2,],ylim=c(0,max(inf)),type='l',col='red') 
polygon(c(rev(dates), dates), 
        c(rev(inf[1,]), inf[3,]), 
        col = rgb(1,0,0,0.3), border = NA)

# plot ppd of npi model fit in CA, TX, FL, NY
# states <- c(5,44,10,35) # CA, TX, FL, NY
load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',2,'-full.Rda'))
log_rt_hat <- sims.npi %>% select(contains('log_rt_hat'))

par(mfrow=c(2,2))
for(state.num in c(5,44,10,35)){
  dates <- ymd(model_data$dates[stan_data$state_inds==state.num])
  abbrev <- model_data$abbrev[state.num]
  load(paste0('state-data/data_',abbrev,'-',state.num,'.Rda'))
  load(paste0('state-sims/sims_',abbrev,'-',state.num,'.Rda'))
  week.inds <- stan_data$error_term_inds[stan.data$state_inds == state.num]
  rt.hat <- apply(
    exp(log_rt_hat[,unique(week.inds)]),2,
    function(x){quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95))}
  )
  rt <- apply(
    sims %>% select(starts_with('r0[')),2,
    function(x){quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95))}
  )
  rt.map <- (sims %>% select(starts_with('r0[')))[which(sims$lp__==max(sims$lp__)),]
  week.dates <- dates[which(c(1,diff(state_data$error_term_inds))==1)]
  
  par(mar=c(4,5,3,1))
  plot(
    week.dates,
    rt.map,
    ylim=range(rt.hat),type='l',col='black',
    ylab=TeX('$R_0(t)$'),xlab = 'Month (2020)',
    main=model_data$states[state.num])  
  abline(h=1,col='grey',lty='dashed')
  polygon(
    c(rev(week.dates), week.dates),
    c(rev(rt.hat[1,]),rt.hat[5,]), 
    col = rgb(0,0,1,0.3), border = NA)
  polygon(
    c(rev(week.dates), week.dates), 
    c(rev(rt.hat[2,]),rt.hat[4,]), 
    col = rgb(0,0,1,0.5), border = NA)
  lines(week.dates,rt.map)  
}





# par(mar=c(4,4,1,0))
plot(dates,rt[3,state_data$error_term_inds],ylim=range(rt),type='l',col='black',
     ylab=TeX('$R_0(t)$'),xlab = 'Month (2020)')  
abline(h=1,col='grey',lty='dashed')
polygon(c(rev(dates), dates), 
        c(rev(rt[1,state_data$error_term_inds]), 
          rt[5,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.3), border = NA)
polygon(c(rev(dates), dates), 
        c(rev(rt[2,state_data$error_term_inds]), 
          rt[4,state_data$error_term_inds]), 
        col = rgb(0,0,1,0.5), border = NA)
lines(dates,rt[3,state_data$error_term_inds])  
