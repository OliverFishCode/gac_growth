# Front-end needs ---------------------------------------------------------
# Load R packages
  library(R2jags)
  library(snowfall)
  library(parallel)
  library(rlecuyer)

  # Need to install JAGS software before you can use this:
  # http://sourceforge.net/projects/mcmc-jags/files/

# Get number of cores on machine
# NOTE that this uses the detectCores() function
# from the 'parallel' package.
  nCpus <- parallel::detectCores() - 1

# Initialize snowfall
  sfInit(parallel = TRUE, cpus=nCpus, type="SOCK")
  
# Wrapper fxn -----
# Define the wrapper function to call in parallel,
# which is really just the whole simulation
  wrapper <- function(x){
    
# . Data definition -----
# Number of age classes
  nages <- 14
  
# Number of samples in each age class
  nsamps <- 30

# Create population of known ages with 30 individuals in each age class
  age <- rep(seq(1, nages, 1), nsamps)
  
# Known parameters of VBGF for simulated population
  linf <- 500
  k <- 0.25
  t0 <- -1.0
  sdlinf <- 0.02
  sdk <- 0.02
  sdt0 <- 0.02
  
  slinf <- rnorm(nages*nsamps, linf, sdlinf)    
  sk <- rnorm(nages*nsamps, k, sdk)             
  st0 <- rnorm(nages*nsamps, t0, sdt0) 
  sw <- slinf*sk  
  
# Add a little random noise to each of the parameters
  slinf <- linf + round(runif(nages*nsamps, -10, 10))
  sk <- k + runif(nages*nsamps, -0.05, 0.05)
  st0 <- t0 + runif(nages*nsamps, -1, 0.5)
  sw <- slinf * sk
  
# Simulated length of individuals based on age and VBGF parameters
  slengq <- (sw/sk)*(1-exp(-sk*(age-st0)))
  
# Put the data together in a dataframe  
  fish <- data.frame(age, slinf, sk, st0, slengq)
  
# . Model definition -----
# Define the model as a function
model <- function(){

  for(i in 1:N){
    # Likelihood
      Y[i] ~ dnorm(L[i], tau[Ti[i]])
      L[i] <- (w[i]/K)*(1-exp(-K*(Ti[i]-to)))
  
    # Linear predictor of w
      log(w[i]) <- beta0
    }
  
  # Priors on VBGF parameters (w defined below by linear model)
    # Brody growth coefficient
      K ~ dunif(0, 1)
    # Age at length zero
      to ~ dunif(-10, 1)
  
  # Priors on parameters of linear model on w
    # Intercept
      beta0 ~ dnorm(0, 0.001)
  
  # Prior distribution for precision at each age
  # This imposes a multiplicative error structure on length at age
    for(t in 1:Tmax){
      tau[t] ~ dgamma(0.01, 0.001)
    }
}
  
  
# . Model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0')
  
# Package the data for JAGS
  vb_data <- list(
    Y = fish$slengq,
    Ti = fish$age,
    Tmax = max(fish$age),
    N = nrow(fish)
  )
 
# Initial values 
  inits <- function(){
    list(
      K = runif(1, 0, 1),
      to = runif(1, -10, 1),
      tau = rgamma(max(fish$age), .01, 1),
      beta0 = rnorm(1, 0, 1)     
    )
  }

# MCMC settings
  ni <- 55000       # Number of draws from posterior (for each chain)
  nt <- 10          # Thinning rate
  nb <- 15000       # Number of draws to discard as burn-in
  nc <- 3           # Number of chains

# Call jags and run the model, re-run if crashes due to
# bad initial values
  vbModgq <- NULL
  attempt <- 0
  while(is.null(vbModgq)){
    attempt <- attempt + 1
    try(
      vbModgq <- jags(data=vb_data, inits=inits, params, model,
        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
        working.directory = getwd(), progress.bar = "none")
    )
  }

# Save the estimates and the known parameter values
  ests <- vbModgq$BUGSoutput$summary[,1] 
  out <- c(ests[-3], linf, k, t0)
  names(out)[4:6] <- c('linf', 'k', 't0')  
  
# Return the result  
  return(list(
    out = out
    ))
  
} # End wrapper function
  
# Load libraries on workers -----
  sfLibrary(R2jags)
  sfLibrary(rlecuyer)
  
# Start network random number generator -----
  sfClusterSetupRNG()
  
# Distribute calculations to workers -----
  # Number of simulations to run
    niterations <- 1000
  
  # Get start time for benchmarking
    start <- Sys.time()
    
  # Run the simulation in parallel
    result <- sfLapply(1:niterations, wrapper)
  
  # Calculate run time
    Sys.time() - start
    
# Stop cluster -----
  sfStop()
  
# Collect results -----
  res <- lapply(result, function(x) x[[c('out')]])
  res <- data.frame(do.call(rbind, res))
    
# Quick sanity check -----
  nrow(res)
  
# Post-processing -----
  par(mar=c(5,5,1,1))
  hist(res$K, col='gray87', xlab=expression(paste(hat(italic('k')))),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(0, .5), main='')
  abline(v=mean(res$k), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=0, las=TRUE)    
