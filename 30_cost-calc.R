rm(list=ls())

# fit linear model for tests administered over time
library(dplyr)
library(readr)
library(lme4)
library(lubridate)

state.pops <- read_csv("data/state_populations_2020.csv")
state.pops <- state.pops %>%
  rename(name = State)

state.abbrev <- read_csv("data/state_abbrev.csv")
state.abbrev <- rbind(state.abbrev,c('District of Columbia','DC'))
state.abbrev <- state.abbrev %>%
  rename(name = state, state = abbreviation)
state.data <- left_join(state.pops,state.abbrev,by='name')

daily.tests <- read_csv("data/daily-tests-per-thousand-people-smoothed-7-day.csv")
daily.tests <- daily.tests %>% 
  filter(Entity == "United States", Day < '2021-01-01') %>%
  rename(tests = new_tests_per_thousand_7day_smoothed)
daily.tests$index <- 1:nrow(daily.tests)
daily.tests$tests <- daily.tests$tests*330000
fit.ow.daily <- lm(tests ~ index, data=daily.tests)
summary(fit.ow.daily)
fit.ow.daily$coefficients[2] # 5822.882 
min(daily.tests$Day) # "2020-03-08"

total.tests <- read_csv("data/full-list-total-tests-for-covid-19.csv")
total.tests <- total.tests %>%
  filter(Entity == "United States", Day < '2021-01-01')
total.tests$index <- 1:nrow(total.tests)
total.tests$index2 <- total.tests$index^2
fit.ow <- lm(total_tests ~ index2, data=total.tests)
summary(fit.ow)
fit.ow$coefficients[2]*2 # 5432.802 
min(total.tests$Day) # "2020-03-01"

total.tests$test.sqrt <- sqrt(total.tests$total_tests)
fit.ow2 <- lm(test.sqrt ~ index, data=total.tests)
summary(fit.ow2)
2*fit.ow2$coefficients[2]^2 # 5731.568 

ctp.data <- read.csv(url("https://covidtracking.com/data/download/all-states-history.csv"))
plot(ctp.data$totalTestResultsIncrease[ctp.data$state=='AK'])
plot(ctp.data$totalTestResults[ctp.data$state=='AK'])
plot(ctp.data$totalTestsViralIncrease[ctp.data$state=='AK'])
plot(ctp.data$totalTestsViral[ctp.data$state=='AK'])
ctp.data <- ctp.data %>%
  filter(date < '2021-01-01', date >= '2020-03-08') %>%
  select(date,state,tests=totalTestResultsIncrease,tot.tests=totalTestResults)
names(ctp.data)
ctp.data$tests.c <- ctp.data$tests - mean(ctp.data$tests)
ctp.data$day <- as.numeric(as.POSIXct(ctp.data$date))
ctp.data$day <- (ctp.data$day - mean(ctp.data$day))/sd(ctp.data$day)
ctp.data$state.fact <- as.factor(ctp.data$state)
fit.lmm <- lmer(tests.c ~ 0+ day + (day|state.fact), data = ctp.data)
summary(fit.lmm)

ctp.data <- left_join(ctp.data,state.data,by='state')
ctp.data <- ctp.data %>%
  filter(!is.na(name))

slopes <- c()
slopes1 <- c()
slopes2 <- c()
r2 <- c()
r2.1 <- c()
r2.2 <- c()
states <- c()
for(state in state.data$state){
  states <- c(states,state)
  y <- 1e6*ctp.data$tests[ctp.data$state == state]/
    state.data$Population[state.data$state == state]
  x <- length(y):1
  slope <- lm(y ~ x, na.action = 'na.omit')
  slopes <- c(slopes, slope$coefficients[2])
  r2 <- c(r2, summary(slope)$r.squared)
  
  y.tot <- 1e6*ctp.data$tot.tests[ctp.data$state == state]/
    state.data$Population[state.data$state == state]
  y.sqrt <- sqrt(y.tot)
  x.sqrt <- length(y.sqrt):1
  x.sq <- x.sqrt^2
  slope <- lm(y.tot ~ x.sq, na.action = 'na.omit')
  slopes1 <- c(slopes1,slope$coefficients[2]*2)
  r2.1 <- c(r2.1, summary(slope)$r.squared)
  
  slope <- lm(y.sqrt ~ x.sqrt, na.action = 'na.omit')
  slopes2 <- c(slopes2,2*slope$coefficients[2]^2)
  r2.2 <- c(r2.2, summary(slope)$r.squared)
}

range(r2.1)
range(r2.2)
mean(r2.1 < r2.2) # slopes2 is a better model
hist(slopes1,breaks=20,xlim=c(0,50)); abline(v=5700/330,col='red')
hist(slopes2,add=TRUE,breaks=20,col=rgb(0,0,1,0.3)); 
state.slopes <- data.frame(state=states,new_tests_per_million_per_day=as.numeric(slopes2))
state.data <- left_join(state.data,state.slopes,by='state')
state.data <- state.data %>% rename(pop2020 = Population)

ctp.usa <- read.csv(url("https://covidtracking.com/data/download/national-history.csv"))
ctp.usa <- ctp.usa %>% 
  filter(date < '2021-01-01', date >= '2020-03-08')
ctp.usa$day <- nrow(ctp.usa):1
rev(ctp.usa$date)[1] # "2020-03-08"


plot(rev(ctp.usa$totalTestResults))
plot(sqrt(rev(ctp.usa$totalTestResults)))
plot(ctp.usa$day^2,ctp.usa$totalTestResults)
plot(rev(ctp.usa$totalTestResultsIncrease))


y <- ctp.usa$totalTestResultsIncrease[ctp.usa$day >= cutoff]
x <- ctp.usa$day[ctp.usa$day >= cutoff]
plot(x,y)
fit.lm <- lm(y ~ x)
summary(fit.lm)
fit.lm$coefficients[2] # 6159.266 

y1 <- sqrt(ctp.usa$totalTestResults[ctp.usa$day >= cutoff])
plot(x,y1)
fit1.lm <- lm(y1 ~ x)
summary(fit1.lm)
2*fit1.lm$coefficients[2]^2 # 5669.259

y2 <- ctp.usa$totalTestResults[ctp.usa$day >= cutoff]
x2 <- x^2
plot(x2,y2)
fit2.lm <- lm(y2 ~ x2) 
fit2.lm$coefficients[2]*2 # 4111.529 

# contact tracing data
# https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2780568

tracing <- data.frame(
  pop = c(
    30781,
    143667,
    623989,
    32149,
    258826,
    293086,
    4217737,
    1934408,
    1110356,
    120629,
    8882190,
    694144,
    884659,
    1138890
  ),
  dept.type = c(
    'local',
    'local',
    'state',
    'local',
    'local',
    'local',
    'state',
    'state',
    'local',
    'local',
    'state',
    'local',
    'state',
    'local'
  ),
  mean.weekly.cases.per100k = c(
    30.9,
    97.1,
    6.3,
    106.5,
    121.2,
    32.2,
    44.3,
    95.1,
    144.5,
    99.3,
    29.3,
    317.6,
    621.9,
    208.6
  ),
  cases = c(
    40,
    589,
    146,
    137,
    718,
    479,
    7041,
    5087,
    7116,
    493,
    10563,
    10757,
    22032,
    8987
  ),
  contacts = c(
    117,
    1146,
    404,
    359,
    712,
    1418,
    10927,
    6068,
    13401,
    173,
    11569,
    2848,
    24190,
    1507
  ),
  prob.interviewed = c(
    1.00,
    0.9915014164305950,
    0.9886685552407930,
    0.9915014164305950,
    0.9121813031161470,
    0.8470254957507080,
    0.8271954674220960,
    0.7790368271954670,
    0.7762039660056660,
    0.7507082152974510,
    0.7393767705382440,
    0.4872521246458920,
    0.4674220963172800,
    0.33144475920679900
  )
)
tracing$interviewed.per100k <- 
  (1e5*tracing$prob.interviewed*tracing$cases/tracing$pop)/4
tracing$mean.weekly.int.per100k <- 
  tracing$mean.weekly.cases.per100k*tracing$prob.interviewed
tracing
tracing %>% filter(dept.type=='state')
sum(tracing$contacts)/sum(tracing$cases) # about 1:1

hist(tracing$mean.weekly.int.per100k,breaks=5)
hist(log(tracing$mean.weekly.int.per100k),breaks=5)

hist(tracing$interviewed.per100k,breaks=5)
hist(log(tracing$interviewed.per100k),breaks=5)
mu <- mean(log(tracing$interviewed.per100k))
sigma <- sd(log(tracing$interviewed.per100k))
y <- (log(tracing$interviewed.per100k) - mu)/sigma
# test for log-normality
shapiro.test(y) # fail to reject
ks.test(y, 'pnorm') # fail to reject

tracing$interviewed.per100k
exp.mean <- exp(mu+sigma^2/2)
exp.sd <- sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
n.tracing <- nrow(tracing)
ci.low <- exp(mu - 1.96*sigma/sqrt(n.tracing)); ci.low #/4
ci.high <- exp(mu + 1.96*sigma/sqrt(n.tracing)); ci.high #/4
exp.mean #/4

# another contact tracing study
# https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2790518
704/(60/7)
895/(60/7)
0.5*895/(60/7)+0.5*704/(60/7)



# state GDP data (BEA)
library(readr)
library(dplyr)
state.gdp <- read_csv("data/SAGDP/SAGDP1__ALL_AREAS_2017_2022.csv")
# Bureau of Economic Analysis 2019 US GDP current dollar growth rate
# https://www.bea.gov/news/2020/gross-domestic-product-fourth-quarter-and-year-2019-third-estimate-corporate-profits
state.gdp <- state.gdp %>%
  filter(Description == 'Current-dollar GDP (millions of current dollars)') %>%
  rename(GDP2019 = `2019`, name = GeoName) %>%
  mutate(GDP2020 = GDP2019*(1.041)) %>%
  select(name,GDP2019,GDP2020)
state.data <- left_join(state.data,state.gdp,by='name')

# state personal income data (BEA)
library(readxl)
state.data$total.income.millions.2019 <- NA
state.data$personal.income2019 <- NA
state.data$pop.2019.bea <- NA
state.data$total.income.millions.2020 <- NA
state.data$personal.income2020 <- NA
state.data$pop.2020.bea <- NA
for(j in 1:nrow(state.data)){
  name <- state.data$name[j]
  state.income <- read_excel("data/covid-workbook-ann.xlsx",sheet = name,skip=4) 
  state.data$total.income.millions.2019[j] <-  state.income[1,3]
  state.data$total.income.millions.2020[j] <- state.income[1,4]
  state.data$pop.2019.bea[j] <- state.income[4,3]
  state.data$pop.2020.bea[j] <- state.income[4,4]
  state.data$personal.income2019[j] <- state.income[5,3]
  state.data$personal.income2020[j] <- state.income[5,4]
}

# personal consumption expenditure / income / GDP / employment data (BEA)
state.consumption <- read_csv("data/SAEXP/SAEXP1__ALL_AREAS_1997_2019.csv")
state.consumption <- state.consumption %>%
  filter(Description == 'Personal consumption expenditures') %>%
  rename(PCE.millions.2019 = `2019`,name = GeoName) %>%
  select(name,PCE.millions.2019)
state.data <- left_join(state.data,state.consumption,by='name')

state.data <- state.data %>%
  rename(GDP.millions.2019 = GDP2019, GDP.millions.2020 = GDP2020)
state.data$PCE2019 <- 1e6*state.data$PCE.millions.2019/as.numeric(state.data$pop.2019.bea)

state.consumption2 <- read_csv("data/SAPCE/SASUMMARY__ALL_AREAS_1998_2022.csv")
state.consumption2 <- state.consumption2 %>%
  rename(name = GeoName)
state.data$GDP.millions.2019.2 <- NA
state.data$GDP.millions.2020.2 <- NA
state.data$total.income.millions.2019.2 <- NA
state.data$total.income.millions.2020.2 <- NA
state.data$PCE.millions.2019.2 <- NA
state.data$PCE.millions.2020.2 <- NA
state.data$personal.income2019.2 <- NA
state.data$personal.income2020.2 <- NA
state.data$PCE2019.2 <- NA
state.data$PCE2020.2 <- NA
state.data$employment2019 <- NA
state.data$employment2020 <- NA

for(j in 1:nrow(state.data)){
  state.name <- state.data$name[j]
  state.consumption.tmp <- state.consumption2 %>%
    filter(name == state.name)
  state.data$GDP.millions.2019.2[j] <- state.consumption.tmp[4,30]
  state.data$GDP.millions.2020.2[j] <- state.consumption.tmp[4,31]
  state.data$total.income.millions.2019.2[j] <- state.consumption.tmp[5,30]
  state.data$total.income.millions.2020.2[j] <- state.consumption.tmp[5,31]
  state.data$PCE.millions.2019.2[j]  <- state.consumption.tmp[7,30]
  state.data$PCE.millions.2020.2[j]  <- state.consumption.tmp[7,31]
  state.data$personal.income2019.2[j]  <- state.consumption.tmp[10,30]
  state.data$personal.income2020.2[j]  <- state.consumption.tmp[10,31]
  state.data$PCE2019.2[j]  <- state.consumption.tmp[12,30]
  state.data$PCE2020.2[j]  <- state.consumption.tmp[12,31] 
  state.data$employment2019[j] <- state.consumption.tmp[15,30]
  state.data$employment2020[j] <- state.consumption.tmp[15,31]
}

state.data <- state.data %>%
  rename(GDP.millions.2020.counterfactual = GDP.millions.2020)
state.data$GDP.percap.2019 <- 
  1e6*as.numeric(state.data$GDP.millions.2019)/
  as.numeric(state.data$pop.2019.bea)
state.data$GDP.percap.2019.2 <- 
  1e6*as.numeric(state.data$GDP.millions.2019.2)/
  as.numeric(state.data$pop.2019.bea)
state.data$GDP.percap.2020.counterfactual <- 
  1e6*as.numeric(state.data$GDP.millions.2020.counterfactual)/
  as.numeric(state.data$pop.2020.bea)
state.data$GDP.percap.2020.2 <- 
  1e6*as.numeric(state.data$GDP.millions.2020.2)/
  as.numeric(state.data$pop.2020.bea)

# other cost parameters
state.data$us.median.income.2019 <- 35980 
state.data$us.median.income.2020 <- 35860
state.data$us.personal.income.2019 <- as.numeric(state.consumption2[10,30][[1]])
state.data$us.personal.income.2020 <- as.numeric(state.consumption2[10,31][[1]])
state.data$median.income.2019.est <- 
  as.numeric(state.data$us.median.income.2019) * 
  as.numeric(state.data$personal.income2019.2)/
  as.numeric(state.data$us.personal.income.2019)
state.data$median.income.2020.est <- 
  as.numeric(state.data$us.median.income.2020) * 
  as.numeric(state.data$personal.income2020.2)/
  as.numeric(state.data$us.personal.income.2020)
state.data$median.income.2020.counterfactual <- 
  state.data$median.income.2019.est * 1.041
state.data$personal.income2020.counterfactual <-
  as.numeric(state.data$personal.income2019.2) * 1.041

# cost parameters
stan_data <- readRDS(file="data/stan_data.RData")
model_data <- readRDS(file="data/model_data.RData")

# daily and weekly discount rates
state.data$daily.discount <- 0.96^(1/366)
state.data$weekly.discount <- 0.96^(1/52)

# cost of fear of infection per capita per week
# Aum et al. (2021) unemployment number
state.data$employment2019.rate <- 
  as.numeric(state.data$employment2019)/as.numeric(state.data$pop.2019.bea)

state.data$cost.fear <-
  stan_data$gamma_mean * # average infection length in days
  0.155 * # approximate case ascertainment rate in SK at this time
  # next two lines are (total employment cost/population)
  state.data$median.income.2020.counterfactual * 0.0268 * 
  1000/366 
# cost per infection per day = (cost/pop)*(pop/inf)/366 = (cost/pop)*1000/366

# VSL, Robinson et al. (2020)
state.data$vsl.low <- 4.47*1e6
state.data$vsl.mid <- 8.31*1e6
state.data$vsl.high <- 10.63*1e6

# productivity cost of an infection: one week worth of wages
state.data$cost.productivity <-
  state.data$median.income.2020.counterfactual/52

# health/medical costs of an infection
# Bartsch et al. (2020)
# and similar to Wharton budget model numbers 
# (https://budgetmodel.wharton.upenn.edu/issues/2020/10/12/covid-trade-offs-in-school-re-opening)
state.data$cost.medical <- 3045

# school closure learning loss costs
# Hanushek and Woessman (2020)
# 69% of GDP per 0.33 years of school closure (4 months)
# 69% of GDP per 0.33 years of school closure 
state.data$cost.learning.high <- 
  3 * 0.69 * state.data$GDP.percap.2020.counterfactual/52

# Psacharopoulos et al. (2021) 
# 9% of GDP per 0.33 years of school closure (4 months)
state.data$cost.learning.low <- 
  3 * 0.09 * state.data$GDP.percap.2020.counterfactual/52

# time until learning loss ends (linearly decreasing cost function)
# assuming a maximum of 0.35 years of learning loss can accrue (Betthäuser et al. (2023))
state.data$time.learning <- round(2*52*0.35) # number of weeks

# direct GDP costs of school closure due to worker absenteeism
# Lempel et al. (2009)
# 0.2% GDP per four weeks of school closure
state.data$cost.schools.direct <- 
  0.002 * state.data$GDP.percap.2020.counterfactual/4

# costs of workplace closure and social distancing
# based on lit review, assume both cause 4 percentage point decline in employment
state.data$cost.work <-
  0.04 * state.data$median.income.2020.counterfactual/52
state.data$cost.social <- state.data$cost.work

# cost of masks
# $0.32 per person per day (see spreadsheet)
# This is one surgical mask per day or one N95 per week per person.
state.data$cost.mask <- 7*0.32

# cost of testing
# $100 per test, 
# multiplied by number of tests added for each week of testing utilization
# assuming testing capacity ramps up at a linear rate when policy in place
state.data$cost.test <-
  100 * 7 * state.data$new_tests_per_million_per_day/1e6
  

# cost of contact tracing
# assume each case costs $66.50 and 
# states ramp up capacity linearly at a rate of 
# 4.72 cases per 100k population per week, based on calculations in spreadsheet
state.data$cost.trace <-
  66.50 * 4.72 / 1e5

for(j in 4:ncol(state.data)){
  print(j)
  state.data[,j] <- as.numeric(state.data[,j][[1]])
}

save(state.data,file='data/state-costs.RDa')
write_csv(state.data,file='data/state-costs.csv')

