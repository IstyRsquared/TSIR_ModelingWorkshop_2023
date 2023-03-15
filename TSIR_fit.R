## Isty  Rysava
## 03/14/2023
## Fitting TSIR simulation model to data - Modeling Methods Workshop

rm(list=ls())

### Data
cases <- read.csv("case_data.csv")

### Sim function - add stochasticity, incursions, vaccination, and seasonality
simTSI_all = function(pop, params, times) {
  
}

### Poisson log-likelihood function
ll_obs_pois <- function(sim, obs){
  L <- sum(dpois(x=obs, lambd=sim+0.0001, log = TRUE))
  return(L)
}

### Estimates
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) 
betas <- seq(20, 30, by=1)
params <- list(beta=NA, alpha=0.97, inc=10/52, vacc = 1/1000, sf=0.3) 
times <- 52*10 
ll <- rep(NA, length(betas)) # initialize empty vector to stor  lll

for(j in 1:length(betas)){
  params$beta <- betas[j]
  out <- simTSI_all(pop, params, times)
  ll[j] <- ll_obs_pois(sim=out$I, obs=cases$cases)
}

### Get best tbeta
betas[which.max(ll)]

### Simulate time series using the fitted beta

### Plot the simulated output over data
