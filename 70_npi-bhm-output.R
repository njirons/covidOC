rm(list=ls())

# aggregate NPI BHM model fits across SEIR posterior trajectories

for(nX in 1:3){ # model with nX confounders
  path <- paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-')
  n.traj <- 100 # number of randomly sampled posterior trajectories from SEIR model
  
  # extract trajectories 
  for(j in 1:n.traj){
    print(j)
    tryCatch({
      load(paste0(path,j,'.Rda'))
      if(j == 1){
        sims.npi <- sims
      }else{
        sims.npi <- rbind(sims.npi,sims)
      }
      rm(sims)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  saveRDS(sims.npi, file=paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-full.rds'))
  save(sims.npi, file=paste0('npi-bhm-sims/npi-bhm-ar1-nX',nX,'-full.Rda'))  
}



