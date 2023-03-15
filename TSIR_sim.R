## Isty  Rysava
## 03/14/2023
## Demonstration of TSIR simulation code - Modeling Methods Workshop

rm(list=ls())

### (1) Deterministic simulation of TSI
## Function
simTSI = function(pop, params, times) {
  lambda = I = S = rep(NA, times) # Initialize empty vectors
  N = pop$N
  B = pop$B
  I[1] = pop$I0 
  S[1] = pop$S0*N
  lambda[1] = pop$I0 
  beta = params$beta
  alpha = params$alpha
  
  for (i in 2:times) {
    lambda[i] = beta * I[i - 1]^alpha * S[i - 1] /N
    if (lambda[i] < 0) {lambda[i] = 0}
    I[i] = lambda[i]
    S[i] = S[i - 1] + B - I[i]
  }
  list(I = I, S = S)
}

## Run
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) # set up initial population list
params <- list(beta = 25, alpha=0.97) # parameters
times <- 52*10 # 10 years, weekly
out <- simTSI(pop, params, times)

## Plot
par(mfrow = c(1, 2))
plot(out$I, ylab = "infected", xlab = "time", type="l")
plot(out$S, out$I, ylab = "infected", xlab = "susceptible", type="l")

### (2) Stochastic simulation of TSI
## Function
simTSI_stoch = function(pop, params, times) {
  lambda = I = S = rep(NA, times) 
  N = pop$N
  B = pop$B
  I[1] = pop$I0 
  S[1] = pop$S0*N
  lambda[1] = pop$I0 
  beta = params$beta
  alpha = params$alpha
  
  for (i in 2:times) {
    lambda[i] = rnorm(1, mean=beta, sd=1) * I[i - 1]^alpha * S[i - 1] /N 
    # consider using Negative Binomial distribution for strong environmental noice
    if (lambda[i] < 0) {lambda[i] = 0}
    I[i] = rpois(1, lambda[i]) 
    S[i] = S[i - 1] + B - I[i]
  }
  list(I = I, S = S)
}

## Run
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) 
params <- list(beta = 25, alpha=0.97) 
times <- 52*10 
out <- simTSI_stoch(pop, params, times)

## Plot
plot(out$I, ylab = "infected", xlab = "time", type="l")
plot(out$S, out$I, ylab = "infected", xlab = "susceptible", type="l")

### (3) Stochastic simulation of TSI with incursions
simTSI_inc = function(pop, params, times) {
  lambda = I = S = rep(NA, times) # Initialize empty vectors
  N = pop$N
  B = pop$B
  I[1] = pop$I0 
  S[1] = pop$S0*N
  lambda[1] = pop$I0 
  beta = params$beta
  alpha = params$alpha
  inc_rate = params$inc
  
  for (i in 2:times) {
    lambda[i] = rnorm(1, mean=beta, sd=1) * I[i - 1]^alpha * S[i - 1] /N
    if (lambda[i] < 0) {lambda[i] = 0}
    I[i] = rpois(1, lambda[i]) + rbinom(1, 1, inc_rate) # incursions can be drawn from Poisson distributions 
    S[i] = S[i - 1] + B - I[i]
  }
  list(I = I, S = S)
}

## Run
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) # set up initial population with 1 infected individual
params <- list(beta = 25, alpha=0.97, inc=10/52) # 10 introductions in a year
times <- 52*10 # 10 years, weekly
out <- simTSI_inc(pop, params, times)

## Plot
plot(out$I, ylab = "infected", xlab = "time", type="l")
plot(out$S, out$I, ylab = "infected", xlab = "susceptible", type="l")

### (4) Stochastic simulation of TSI with lasting immunity due to vaccination
simTSI_vacc = function(pop, params, times) {
  lambda = I = S = rep(NA, times) # Initialize empty vectors
  N = pop$N
  B = pop$B
  I[1] = pop$I0 
  S[1] = pop$S0*N
  lambda[1] = pop$I0 
  beta = params$beta
  alpha = params$alpha
  vacc_rate = params$vacc
  
  for (i in 2:times) {
    lambda[i] = rnorm(1, mean=beta, sd=1) * I[i - 1]^alpha * S[i - 1] /N
    if (lambda[i] < 0) {lambda[i] = 0}
    I[i] = rpois(1, lambda[i]) 
    S[i] = S[i - 1] + B - I[i] -  rbinom(1, S[i - 1], vacc_rate)
  }
  list(I = I, S = S)
}

## Run
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) 
params <- list(beta = 25, alpha=0.97, vacc = 1/1000) 
times <- 52*10 
out <- simTSI_vacc(pop, params, times)

## Plot
plot(out$I, ylab = "infected", xlab = "time", type="l")
plot(out$S, out$I, ylab = "infected", xlab = "susceptible", type="l")

### Stochastic simulation of TSI with seasonality
simTSI_seas = function(pop, params, times) {
  lambda = I = S = rep(NA, times) # Initialize empty vectors
  N = pop$N
  B = pop$B
  I[1] = pop$I0 
  S[1] = pop$S0*N
  lambda[1] = pop$I0 
  sf  = params$sf
  beta = params$beta*(1 + sf*cos(2*pi*(1:(times))/52))
  alpha = params$alpha
  
  for (i in 2:times) {
    lambda[i] = rnorm(1, mean=beta[i], sd=1) * I[i - 1]^alpha * S[i - 1] /N
    if (lambda[i] < 0) {lambda[i] = 0}
    I[i] = rpois(1, lambda[i])  
    S[i] = S[i - 1] + B - I[i] 
  }
  list(I = I, S = S)
}

## Run
pop <- list(S0=0.06, N=3.3E6, I0=180, B=2300) 
params <- list(beta = 25, alpha=0.97, sf=0.2) 
times <- 52*10 
out <- simTSI_seas(pop, params, times)

## Plot
plot(out$I, ylab = "infected", xlab = "time", type="l")
plot(out$S, out$I, ylab = "infected", xlab = "susceptible", type="l")

