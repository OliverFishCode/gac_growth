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
  sfInit(parallel = TRUE, cpus=ncpus, type="SOCK")

# Wrapper fxn -----
# Define the wrapper function to call in parallel,
# which is really just the whole simulation
  wrapper <- function(x){
    
# . Environment variables
  #Sys.setenv(PATH=$PATH:/util/academic/jags/4.2.0/bin/jags)
    
# . Data definition -----
# Number of age classes
  nages = 10
  
# Number of samples in each age class
  nsamps = 10

# Create population of known ages with nsamps
# individuals in each age class
  age = rep(seq(1, nages, 1), nsamps)
  
# Define parameters
  k = 0.3
  w = 100
  t0 = -1  
  
# Randomly sample parameters for each fish
# K is on the log scale for simulation to 
# constrain real value of samples to 
# be positive.
  sampl <- 
  rnorm(n=nages*nsamps,
        mean=(w/k)*(1-exp(-k*(age-t0))),
        sd=10)
  
# Combine the von Bert parameters from each
# sample with the age vector
  fish <- data.frame(age, k, w, t0, sampl)
  names(fish)[2:5] <- c("sk", "sw", "st0", "length")
  
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
        to ~ dunif(-10, 10)
    # Priors on parameters of linear model on w
      # Intercept
        beta0 ~ dnorm(0, 0.001)
    # Prior distribution for precision at each age
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
    Y = fish$length,
    Ti = fish$age,
    Tmax = max(fish$age),
    N = nrow(fish)
  )
 
# Initial values 
  inits <- function(){
    list(
      K = runif(1, 0, 1),
      to = runif(1, -10, 10),
      tau = rgamma(max(fish$age), .01, 1),
      beta0 = rnorm(1, 0, 1)     
    )
  }

# MCMC settings
  ni <- 55000        # Number of samples
  nt <- 10           # Thinning rate
  nb <- 15000        # Burn-in
  nc <- 3            # Number of chains

# Call jags and run the model, re-run if crashes
# due to bad initial values
  fin <- NULL
  attempt <- 0
  while(is.null(fin)){
    attempt <- attempt + 1
    try({
      vbMod <- jags.model(file=textConnection(modelString),
                          data=vb_data,
                          inits=inits,
                          n.adapt = 100,
                          quiet=TRUE)
      update(vbMod, nb, progress.bar = "none")
      fin <- coda.samples(vbMod, params, ni, thin=nt,
                          progress.bar = "none")
      }
    )
  }

# Save the estimates and the known parameter values
  ests <- summary(fin)[[1]][,1]
  out <- unlist(c(ests, fish[1,2:4]))

# Return the result  
  return(list(
    out = out
    ))
  
} # End wrapper function
  
# Load libraries on workers -----
  sfLibrary(rjags)
  sfLibrary(rlecuyer)
  
# Start network random number generator -----
  #sfClusterSetupRNG()
  
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
    
# Save results
  save(res, file="result.rda")
  