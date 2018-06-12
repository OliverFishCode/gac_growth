# Front-end needs ---------------------------------------------------------
# Load R packages
  library(rjags)
  library(snowfall)
  library(rlecuyer)

# Get number of cores
args = commandArgs(trailingOnly = TRUE);
ncpus = args[1];
#ncpus = 3 # Uncomment to run on local workstation

# Initialize snowfall
<<<<<<< HEAD:parallelTesting/basemodelParallel.R
  sfInit(parallel = TRUE, cpus=nCpus, type="SOCK")
  
=======
  sfInit(parallel = TRUE, cpus=ncpus, type="SOCK")

>>>>>>> 69a0bdc76e9e64c7bf38cf4b33c8ce7d3b42ff18:parallelTesting/gacbasemodeltest.R
# Wrapper fxn -----
# Define the wrapper function to call in parallel,
# which is really just the whole simulation
  wrapper <- function(x){
    
# . Environment variables
  # Sys.setenv(PATH=$PATH:/usr/lib/rjags.so)
    
# . Data definition -----
# Number of age classes
  nages <- 10
  
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
  modelString ="
    model{
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
  "
  
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
  ni <- 550       # Number of draws from posterior (for each chain)
  nt <- 10          # Thinning rate
  nb <- 150       # Number of draws to discard as burn-in
  nc <- 3           # Number of chains

# Call jags and run the model, re-run if crashes due to
# bad initial values
  fin <- NULL
  attempt <- 0
  while(is.null(fin)){
    attempt <- attempt + 1
    try({
      vbMod <- jags.model(file=textConnection(modelString) , data=vb_data,
                            inits=inits, n.adapt = 100, quiet=TRUE)
      update(vbMod, nb, progress.bar = "none")
      fin <- coda.samples(vbMod, params, ni, progress.bar = "none")
      }

    )
  }

# Save the estimates and the known parameter values
  ests <- summary(fin)[[1]][,1]
  out <- c(ests, linf, k, t0)
  names(out)[4:6] <- c('linf', 'k', 't0')  
  
# Return the result  
  return(list(
    out = out
    ))
  
} # End wrapper function
  
# Load libraries on workers -----
  sfLibrary(rjags)
  sfLibrary(rlecuyer)
  
# Start network random number generator -----
  sfClusterSetupRNG()
  
# Distribute calculations to workers -----
  # Number of simulations to run
    niterations <- 48
  
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
    
# Save results
  save(res, file="result.rda")
  