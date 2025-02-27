# calculate rt based on parameter estimates and NPI policy
rt.eval <- function(u,u_mean,X,old.residue,eps,mu0,mu_npi,mu_X,phi.ar){
  if(is.null(ncol(mu_X))){
    lin.pred <- mu0 + mu_npi %*% (u - u_mean) + mu_X * X
  }else{
    lin.pred <- mu0 + mu_npi %*% (u - u_mean) + colSums(t(mu_X) * X) 
  }
  log_rt <- lin.pred + phi.ar * old.residue + eps
  residue <- log_rt - lin.pred
  
  return(
    list(rt = exp(log_rt), residue = residue)
  )
}

# simulate SEIRD trajectories under NPI policy
seird.sim <- function(u, params){
  n.days <- params$n.days
  n.weeks <- params$n.weeks
  week.ind <- params$week.ind
  n.iter <- params$n.iter
  ifr <- params$ifr
  gamma <- params$gamma
  delta <- params$delta
  mu <- params$mu
  rr <- params$rr
  n.pop <- params$n.pop
  d0 <- params$d0
  mu0 <- params$mu0
  mu_npi <- params$mu_npi
  mu_X <- params$mu_X
  phi.ar <- params$phi.ar
  eps <- params$eps
  d_mult <- params$d_mult
  c_mult <- params$c_mult
  u_mean <- params$u_mean
  nX <- params$nX
  
  sus <- matrix(NA,nrow=n.iter,ncol=n.days)
  exps <- matrix(NA,nrow=n.iter,ncol=n.days)
  inf <- matrix(NA,nrow=n.iter,ncol=n.days)
  rem_s <- matrix(NA,nrow=n.iter,ncol=n.days)
  rem_d <- matrix(NA,nrow=n.iter,ncol=n.days)
  d <- matrix(NA,nrow=n.iter,ncol=n.days)
  exps.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
  inf.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
  rem_d.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
  rt.week <- matrix(NA,nrow=n.iter,ncol=n.weeks)  
  residue.week <- matrix(NA,nrow=n.iter,ncol=n.weeks)
  nu <- matrix(NA,nrow=n.iter,ncol=n.days)
  
  # define SEIRD states
  sus[,1] <- params$sus0
  exps[,1] <- params$exps0
  inf[,1] <- params$inf0
  rem_s[,1] <- params$rem_s0
  rem_d[,1] <- params$rem_d0
  d[,1] <- pmin(1,pmax(0,rem_d[,1]*mu*n.pop))
  exps.week0 <- params$inf0
  inf.week0 <- params$rem_d0
  rem_d.week0 <- params$d0
  exps.week [,week.ind[1]] <- exps.week[,week.ind[1]] + exps[,1]
  inf.week [,week.ind[1]] <- inf.week[,week.ind[1]] + inf[,1]
  rem_d.week[,week.ind[1]] <- rem_d.week[,week.ind[1]] + rem_d[,1]
  X <-rbind(
    d_mult * rem_d.week0,
    c_mult * inf.week0,
    c_mult * exps.week0
  )
  
  # calculate new infections
  rt.sim <- rt.eval(
    u=u[,week.ind[1]],
    u_mean = u_mean,
    X = X[1:nX,],
    old.residue=0,
    eps=eps[,week.ind[1]],
    mu0 = mu0,
    mu_npi = mu_npi,
    mu_X = mu_X,
    phi.ar = phi.ar
  )
  rt.week[,week.ind[1]] <- rt.sim$rt
  residue.week[,week.ind[1]] <- rt.sim$residue
  nu[,1] <- pmin(1,pmax(0,gamma*rt.week[,week.ind[1]]*sus[,1]*inf[,1]))  
  
  for(j in 1:(n.days-1)){
    # evolve SEIRD states forward in time
    sus[,j+1] <- pmin(1,pmax(0,sus[,j] + rr*rem_s[,j] - nu[,j]))
    exps[,j+1] <- pmin(1,pmax(0,exps[,j]*(1-delta) + nu[,j]))
    inf[,j+1] <- pmin(1,pmax(0,inf[,j]*(1-gamma) + delta*exps[,j]))
    rem_s[,j+1] <- pmin(1,pmax(0,rem_s[,j]*(1-rr) + gamma*(1-ifr)*inf[,j]))
    rem_d[,j+1] <- pmin(1,pmax(0,rem_d[,j]*(1-mu) + gamma*ifr*inf[,j]))
    d[,j+1] <- n.pop * pmin(1,pmax(0,rem_d[,j+1]*mu))  
    exps.week [,week.ind[j+1]] <- exps.week[,week.ind[j+1]] + exps[,j+1] 
    inf.week [,week.ind[j+1]] <- inf.week[,week.ind[j+1]] + inf[,j+1] 
    rem_d.week [,week.ind[j+1]] <- rem_d.week[,week.ind[j+1]] + rem_d[,j+1]
    
    # calculate new infections
    if(week.ind[j+1] > week.ind[j]){
      X <-rbind(
        d_mult * rem_d.week[,week.ind[j+1]-1],
        c_mult * inf.week[,week.ind[j+1]-1],
        c_mult * exps.week[,week.ind[j+1]-1]
      )
      
      rt.sim <- rt.eval(
        u=u[,week.ind[j+1]],
        u_mean = u_mean,        
        X = X[1:nX,],
        old.residue=residue.week[,week.ind[j+1]-1],
        eps=eps[,week.ind[j+1]],
        mu0 = mu0,
        mu_npi = mu_npi,
        mu_X = mu_X,
        phi.ar = phi.ar
      )
      rt.week[,week.ind[j+1]] <- rt.sim$rt
      residue.week[,week.ind[j+1]] <- rt.sim$residue      
    }
    nu[,j+1] <- pmin(1,pmax(0,gamma*rt.week[,week.ind[j+1]]*sus[,j+1]*inf[,j+1]))
  }
  return(list(
    sus = sus, exps = exps, inf = inf, rem_s = rem_s, rem_d = rem_d, nu = nu, d = d, 
    exps.week = exps.week, inf.week = inf.week, rem_d.week = rem_d.week, 
    rt.week = rt.week, residue.week=residue.week
  ))
}

# define cost function
cost.schools.opt <- function(sum.school,params){
  tot.learning <- min(sum.school/params$time.learning,1)
  return(
    sum.school * params$cost.schools.direct +
      0.5 * params$cost.learning * params$time.learning * (1 - (1-tot.learning)^2)
  )
}
cost.schools <- function(sum.school,params){
  tot.school <- max(0,sum.school - sum(params$normal.school.year)) # number of weeks w/o school in a normal year
  tot.learning <- min(tot.school/params$time.learning,1)
  return(
    tot.school * params$cost.schools.direct +
      0.5 * params$cost.learning * params$time.learning * (1 - (1-tot.learning)^2)
  )
}

cost.test <- function(u.test,u.test.past,params){
  (u.test.past+u.test)*params$cost.test
}
cost.trace <- function(u.trace,u.trace.past,params){
  (u.trace.past+u.trace)*params$cost.trace
}

make.npis <- function(u,params){
  p <- nrow(u)
  if(p == 3){
    u <- rbind(
      u[1:3,],u[3,],u[3,],u[3,],u[3,],u[3,],
      matrix(1,nrow=3,ncol=params$n.weeks))
  }else if(p == 6){
    u <- rbind(u[1:3,],u[3,],u[3,],u[3,],u[3,],u[3:6,])
  }else if(p != 11){
    return('Error: NPI matrix is of incorrect dimension (i.e., not 3, 6, or 11).')
  }
  return(u)
}

cost.opt <- function(u.vec,params){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools.opt(sum(u[1,]),params)
  cost.work.total <- params$cost.work * sum(u[2,])
  cost.social.total <- params$cost.social * sum(colMeans(u[3:8,]))
  cost.test.total <- sum(
    cost.test(u[9,],cumsum(c(0,u[9,1:(params$n.weeks-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,],cumsum(c(0,u[10,1:(params$n.weeks-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * sum(colMeans(traj.sim$nu))
  
  cost.val <- cost.nu + 
    cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  return(cost.val)
}

cost.final <- function(u.vec,params){
  u.school <- params$normal.school.year
  u.school[u.school == 0] <- u.vec[1:params$n.school.weeks]
  u <- make.npis(matrix(
    c(u.school,u.vec[-(1:params$n.school.weeks)]),byrow=TRUE,
    ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,]),params)
  cost.work.total <- params$cost.work * sum(u[2,])
  cost.social.total <- params$cost.social * sum(colMeans(u[3:8,]))
  cost.test.total <- sum(
    cost.test(u[9,],cumsum(c(0,u[9,1:(params$n.weeks-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,],cumsum(c(0,u[10,1:(params$n.weeks-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * sum(colMeans(traj.sim$nu))
  
  cost.val <- cost.nu + 
    cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  return(cost.val)
}

cost.workplace <- function(u.vec,params){
  u.school <- params$normal.school.year
  u <- make.npis(matrix(
    c(u.school,u.vec,rep(1,4*params$n.weeks)),byrow=TRUE,
    ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,]),params)
  cost.work.total <- params$cost.work * sum(u[2,])
  cost.social.total <- params$cost.social * sum(colMeans(u[3:8,]))
  cost.test.total <- sum(
    cost.test(u[9,],cumsum(c(0,u[9,1:(params$n.weeks-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,],cumsum(c(0,u[10,1:(params$n.weeks-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * sum(colMeans(traj.sim$nu))
  
  cost.val <- cost.nu + 
    cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  return(cost.val)
}

cost <- function(u.vec,params){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,]),params)
  cost.work.total <- params$cost.work * sum(u[2,])
  cost.social.total <- params$cost.social * sum(colMeans(u[3:8,]))
  cost.test.total <- sum(
    cost.test(u[9,],cumsum(c(0,u[9,1:(params$n.weeks-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,],cumsum(c(0,u[10,1:(params$n.weeks-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * sum(colMeans(traj.sim$nu))
  
  cost.val <- cost.nu + 
    cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  return(cost.val)
}

cost.sim <- function(u.vec,nu.sims,params){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,]),params)
  cost.work.total <- params$cost.work * sum(u[2,])
  cost.social.total <- params$cost.social * sum(colMeans(u[3:8,]))
  cost.test.total <- sum(
    cost.test(u[9,],cumsum(c(0,u[9,1:(params$n.weeks-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,],cumsum(c(0,u[10,1:(params$n.weeks-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,]) 
  
  cost.nu <- params$cost.nu * sum(colMeans(nu.sims))
  
  cost.val <- cost.nu + 
    cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  return(cost.val)
}

cost.week <- function(u.vec,params,week){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,1:week]),params)
  cost.work.total <- params$cost.work * sum(u[2,week])
  cost.social.total <- params$cost.social * sum(u[3:8,week])
  cost.test.total <- sum(
    cost.test(u[9,week],cumsum(c(0,u[9,1:(week-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,week],cumsum(c(0,u[10,1:(week-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,week]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * 
    sum(colMeans(traj.sim$nu[,which(params$week.ind == week)]))
  
  cost.npis <- cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  
  return(list(cost.npis = cost.npis, cost.inf = cost.nu))
}

cost.cumulative <- function(u.vec,params,week){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,1:week]),params)
  cost.work.total <- params$cost.work * sum(u[2,1:week])
  if(week > 1){
    cost.social.total <- params$cost.social * sum(colMeans(u[3:8,1:week])) 
  }else{
    cost.social.total <- params$cost.social * mean(u[3:8,week]) 
  }
  cost.test.total <- sum(
    cost.test(u[9,1:week],cumsum(c(0,u[9,1:(week-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,1:week],cumsum(c(0,u[10,1:(week-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,1:week]) 
  
  traj.sim <- seird.sim(u,params)
  cost.nu <- params$cost.nu * 
    rowSums(traj.sim$nu[,which(params$week.ind <= week)])
  
  cost.npis <- cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  
  return(list(cost.npis = cost.npis, cost.inf = cost.nu))
}

cost.week.sim <- function(u.vec,nu.sims,params,week){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,1:week]),params)
  cost.work.total <- params$cost.work * sum(u[2,week])
  cost.social.total <- params$cost.social * sum(u[3:8,week])
  cost.test.total <- sum(
    cost.test(u[9,week],cumsum(c(0,u[9,1:(week-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,week],cumsum(c(0,u[10,1:(week-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,week]) 
  
  cost.nu <- params$cost.nu * 
    sum(colMeans(nu.sims[,which(params$week.ind == week)]))
  
  cost.npis <- cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  
  return(list(cost.npis = cost.npis, cost.inf = cost.nu))
}

cost.cumulative.sim <- function(u.vec,nu.sims,params,week){
  u <- make.npis(matrix(u.vec,ncol = params$n.weeks),params)
  
  cost.schools.total <- cost.schools(sum(u[1,1:week]),params)
  cost.work.total <- params$cost.work * sum(u[2,1:week])
  if(week > 1){
    cost.social.total <- params$cost.social * sum(colMeans(u[3:8,1:week])) 
  }else{
    cost.social.total <- params$cost.social * mean(u[3:8,week]) 
  }
  cost.test.total <- sum(
    cost.test(u[9,1:week],cumsum(c(0,u[9,1:(week-1)])),params)
  )
  cost.trace.total <- sum(
    cost.trace(u[10,1:week],cumsum(c(0,u[10,1:(week-1)])),params)
  )  
  cost.mask.total <- params$cost.mask * sum(u[11,1:week]) 
  
  cost.nu <- params$cost.nu * 
    rowSums(nu.sims[,which(params$week.ind <= week)])
  
  cost.npis <- cost.schools.total +
    cost.work.total + 
    cost.social.total +
    cost.test.total +
    cost.trace.total +
    cost.mask.total
  
  return(list(cost.npis = cost.npis, cost.inf = cost.nu))
}

pacman::p_load(optimParallel,snow)
nCores <- detectCores(logical = FALSE); nCores
nThreads<- detectCores(logical = TRUE); nThreads
cluster <- makeCluster(4, type = "SOCK")
setDefaultCluster(cluster)
clusterEvalQ(
  cluster, 
  make.npis <- function(u,params){
    p <- nrow(u)
    if(p == 3){
      u <- rbind(
        u[1:3,],u[3,],u[3,],u[3,],u[3,],u[3,],
        matrix(1,nrow=3,ncol=params$n.weeks))
    }else if(p == 6){
      u <- rbind(u[1:3,],u[3,],u[3,],u[3,],u[3,],u[3:6,])
    }else if(p != 11){
      return('Error: NPI matrix is of incorrect dimension (i.e., not 3, 6, or 11).')
    }
    return(u)
  }
)
clusterEvalQ(
  cluster, 
  cost.schools <- function(sum.school,params){
    tot.school <- max(0,sum.school - sum(params$normal.school.year)) # number of weeks w/o school in a normal year
    tot.learning <- min(tot.school/params$time.learning,1)
    return(
      tot.school * params$cost.schools.direct +
        0.5 * params$cost.learning * params$time.learning * (1 - (1-tot.learning)^2)
    )
  }
)
clusterEvalQ(
  cluster, 
  cost.test <- function(u.test,u.test.past,params){
    (u.test.past+u.test)*params$cost.test
  }
)
clusterEvalQ(
  cluster, 
  cost.trace <- function(u.trace,u.trace.past,params){
    (u.trace.past+u.trace)*params$cost.trace
  }
)
clusterEvalQ(
  cluster, 
  rt.eval <- function(u,u_mean,X,old.residue,eps,mu0,mu_npi,mu_X,phi.ar){
    if(is.null(ncol(mu_X))){
      lin.pred <- mu0 + mu_npi %*% (u - u_mean) + mu_X * X    
    }else{
      lin.pred <- mu0 + mu_npi %*% (u - u_mean) + colSums(t(mu_X) * X) 
    }
    log_rt <- lin.pred + phi.ar * old.residue + eps
    residue <- log_rt - lin.pred
    
    return(
      list(rt = exp(log_rt), residue = residue)
    )
  }
)
clusterEvalQ(
  cluster,
  seird.sim <- function(u, params){
    n.days <- params$n.days
    n.weeks <- params$n.weeks
    week.ind <- params$week.ind
    n.iter <- params$n.iter
    ifr <- params$ifr
    gamma <- params$gamma
    delta <- params$delta
    mu <- params$mu
    rr <- params$rr
    n.pop <- params$n.pop
    d0 <- params$d0
    mu0 <- params$mu0
    mu_npi <- params$mu_npi
    mu_X <- params$mu_X
    phi.ar <- params$phi.ar
    eps <- params$eps
    d_mult <- params$d_mult
    c_mult <- params$c_mult
    u_mean <- params$u_mean
    nX <- params$nX
    
    sus <- matrix(NA,nrow=n.iter,ncol=n.days)
    exps <- matrix(NA,nrow=n.iter,ncol=n.days)
    inf <- matrix(NA,nrow=n.iter,ncol=n.days)
    rem_s <- matrix(NA,nrow=n.iter,ncol=n.days)
    rem_d <- matrix(NA,nrow=n.iter,ncol=n.days)
    d <- matrix(NA,nrow=n.iter,ncol=n.days)
    exps.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
    inf.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
    rem_d.week <- matrix(0,nrow=n.iter,ncol=n.weeks)
    rt.week <- matrix(NA,nrow=n.iter,ncol=n.weeks)  
    residue.week <- matrix(NA,nrow=n.iter,ncol=n.weeks)
    nu <- matrix(NA,nrow=n.iter,ncol=n.days)
    
    # define SEIRD states
    sus[,1] <- params$sus0
    exps[,1] <- params$exps0
    inf[,1] <- params$inf0
    rem_s[,1] <- params$rem_s0
    rem_d[,1] <- params$rem_d0
    d[,1] <- pmin(1,pmax(0,rem_d[,1]*mu*n.pop))
    exps.week0 <- params$inf0
    inf.week0 <- params$rem_d0
    rem_d.week0 <- params$d0
    exps.week [,week.ind[1]] <- exps.week[,week.ind[1]] + exps[,1]
    inf.week [,week.ind[1]] <- inf.week[,week.ind[1]] + inf[,1]
    rem_d.week[,week.ind[1]] <- rem_d.week[,week.ind[1]] + rem_d[,1]
    X <-rbind(
      d_mult * rem_d.week0,
      c_mult * inf.week0,
      c_mult * exps.week0
    )
    
    # calculate new infections
    rt.sim <- rt.eval(
      u=u[,week.ind[1]],
      u_mean = u_mean,
      X = X[1:nX,],
      old.residue=0,
      eps=eps[,week.ind[1]],
      mu0 = mu0,
      mu_npi = mu_npi,
      mu_X = mu_X,
      phi.ar = phi.ar
    )
    rt.week[,week.ind[1]] <- rt.sim$rt
    residue.week[,week.ind[1]] <- rt.sim$residue
    nu[,1] <- pmin(1,pmax(0,gamma*rt.week[,week.ind[1]]*sus[,1]*inf[,1]))  
    
    for(j in 1:(n.days-1)){
      # evolve SEIRD states forward in time
      sus[,j+1] <- pmin(1,pmax(0,sus[,j] + rr*rem_s[,j] - nu[,j]))
      exps[,j+1] <- pmin(1,pmax(0,exps[,j]*(1-delta) + nu[,j]))
      inf[,j+1] <- pmin(1,pmax(0,inf[,j]*(1-gamma) + delta*exps[,j]))
      rem_s[,j+1] <- pmin(1,pmax(0,rem_s[,j]*(1-rr) + gamma*(1-ifr)*inf[,j]))
      rem_d[,j+1] <- pmin(1,pmax(0,rem_d[,j]*(1-mu) + gamma*ifr*inf[,j]))
      d[,j+1] <- n.pop * pmin(1,pmax(0,rem_d[,j+1]*mu))  
      exps.week [,week.ind[j+1]] <- exps.week[,week.ind[j+1]] + exps[,j+1] 
      inf.week [,week.ind[j+1]] <- inf.week[,week.ind[j+1]] + inf[,j+1] 
      rem_d.week [,week.ind[j+1]] <- rem_d.week[,week.ind[j+1]] + rem_d[,j+1]
      
      # calculate new infections
      if(week.ind[j+1] > week.ind[j]){
        X <-rbind(
          d_mult * rem_d.week[,week.ind[j+1]-1],
          c_mult * inf.week[,week.ind[j+1]-1],
          c_mult * exps.week[,week.ind[j+1]-1]
        )
        
        rt.sim <- rt.eval(
          u=u[,week.ind[j+1]],
          u_mean = u_mean,        
          X = X[1:nX,],
          old.residue=residue.week[,week.ind[j+1]-1],
          eps=eps[,week.ind[j+1]],
          mu0 = mu0,
          mu_npi = mu_npi,
          mu_X = mu_X,
          phi.ar = phi.ar
        )
        rt.week[,week.ind[j+1]] <- rt.sim$rt
        residue.week[,week.ind[j+1]] <- rt.sim$residue      
      }
      nu[,j+1] <- pmin(1,pmax(0,gamma*rt.week[,week.ind[j+1]]*sus[,j+1]*inf[,j+1]))
    }
    return(list(
      sus = sus, exps = exps, inf = inf, rem_s = rem_s, rem_d = rem_d, nu = nu, d = d, 
      exps.week = exps.week, inf.week = inf.week, rem_d.week = rem_d.week, 
      rt.week = rt.week, residue.week=residue.week
    ))
  } 
)