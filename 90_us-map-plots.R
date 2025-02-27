rm(list=ls())

# plot US state map with colors determined by total NPI effect / fear of infections
pacman::p_load(
  ggplot2,rjson,readr,dplyr,tidyr,cmdstanr,maps,mapproj,latex2exp
  )
require(maps)

states_map <- map_data("state")

# read in data
stan.data <- readRDS(file="data/stan_data.RData")
model.data <- readRDS(file="data/model_data.RData")

# number of NPIs
p <- 11

# load NPI regression model estimates
load(paste0('npi-bhm-sims/npi-bhm-ar1-nX',2,'-full.Rda'))

r0 <- rep(NA,stan.data$s)
total.effect <- rep(NA,stan.data$s)
infection.fear <- rep(NA,stan.data$s)
for(state in 1:stan.data$s){
  print(state)
  
  r0[state] <- median((sims.npi %>% select(paste0('R0_s[',state,']')))[[1]])
  
  mu_s <- as.matrix(sims.npi %>% select(starts_with(paste0('mu_s[',state,','))))
  mu_npi <- mu_s[,2:(p+1)]
  total.effect[state] <- 1-exp(median(rowSums(mu_npi)))
  
  infection.fear[state] <- 1-exp(median(mu_s[,p+2]))
  
  rm(mu_s,mu_npi)
}

df <- data.frame(
  region = tolower(model.data$states),
  R0 = r0, `NPI effect` = 100*total.effect, 
  `Death effect` = 100*infection.fear
  )

df.map <- left_join(states_map, df, by = "region")

map.plot <- ggplot()+
  geom_polygon(
    data=df.map,
    aes(x=long, y=lat, group = group, fill=R0),color="grey50")+
  labs(x="",y="",title=TeX("Basic reproduction number $R_0$"))+
  theme_void() +  
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  scale_fill_gradient(low="#B9BDF1", high="#4679A3")
map.plot

map.plot <- ggplot()+
  geom_polygon(
    data=df.map,
    aes(x=long, y=lat, group = group, fill=NPI.effect),color="grey50")+
  labs(x="",y="",title=TeX("Percent reduction in $R_0$ due to NPIs"))+
  theme_void() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
    ) +
  scale_fill_gradient(low="#B9BDF1", high="#4679A3")
map.plot

map.plot <- ggplot()+
  geom_polygon(
    data=df.map,
    aes(x=long, y=lat, group = group, fill=Death.effect),color="grey50")+
  labs(x="",y="",title=TeX("Percent reduction in $R_0$ per COVID death per 10,000 population"))+
  theme_void() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  scale_fill_gradient(low="#B9BDF1", high="#4679A3")
map.plot



