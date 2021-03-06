# Front-end needs ---------------------------------------------------------
# Load R packages
  library(rjags)
  library(snowfall)
  library(rlecuyer)

# Get number of cores
args = commandArgs(trailingOnly = TRUE);
ncpus = args[1];
#ncpus = 3  # Uncomment to run on local workstation

# Initialize snowfall
  sfInit(parallel = TRUE, cpus=ncpus, type="SOCK")

# Wrapper fxn -----
# Define the wrapper function to call in parallel,
# which is really just the whole simulation
  wrapper <- function(x){
    
# . Data definition -----
# Number of age classes
  nages = 10
  
# Number of samples in each age class
  nsamps = 20

# Create population of known ages with 30 individuals in each age class
  age = rep(seq(1, nages, 1), nsamps)
  
# Known parameters of VBGF for simulated population
# Define parameters
  sk = 0.3
  st0 = -1  
  sbeta0 = 5     # Mean of intercept for glm on omega
  sbetaT = .12   # Mean effect of temperature
  temp = runif(nages*nsamps, 15, 25) # Simulated temperatures
  sw = exp(sbeta0 + sbetaT*as.vector(scale(temp)))   
  
# Simulated length of individuals based on age and VBGF parameters
  slengq = rnorm(nages*nsamps,
                 (sw/sk)*(1-exp(-sk*(age-st0))),
                 10)

# Put the simulated data together in a dataframe
  fish = data.frame(age,    # Fish age
                    temp,   # Temperature
                    sw,     # Omega derived
                    sbeta0, # Omega intercept
                    sbetaT, # Omega slope
                    sk,     # Brody growth coefficient
                    st0,    # t0
                    slengq  # Length from gq params and age- identical
                    )
  
# . Model definition ----- 
# Define the model as a function
modelString = "
    model{
      # Model and likelihood definitions, row based
        for(i in 1:N){
        # Likelihood
          Y[i] ~ dnorm(L[i], tau[Ti[i]])

        # von Bertalanffy growth function
          L[i] <- (w[i]/K)*(1-exp(-K*(Ti[i]-to)))

        # Linear predictor of omega
          log(w[i]) <- beta0 + betaT*temp[i]
        }

      # Priors on VBGF parameters 
      # omega defined above by linear model
        # Brody growth coefficient
          K ~ dunif(0, 1)
        # Age at length zero
          to ~ dunif(-10, 1)

      # Priors on parameters of general linearized model on w
        # Intercept
          beta0 ~ dnorm(0, 0.001)
        # Effect of temperature
          betaT ~ dnorm(0, 0.001)

      # Prior distribution for precision at each age
      # This imposes a multiplicative error structure on length at age
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.01, 0.001)
        }
    }"
  
# . Model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0', 'betaT')
  
# Package the data for JAGS
  vb_data = list(
    Y = fish$slengq,
    Ti = fish$age,
    Tmax = max(fish$age),
    temp = as.vector(scale(fish$temp)),
    N = nrow(fish)
  )
 
# Initial values 
  inits <- function(){
    list(
      K = runif(1, 0, 1),
      to = runif(1, -10, 1),
      tau = rgamma(max(fish$age), .01, 1),
      beta0 = rnorm(1, 0, 1),
      betaT = rnorm(1, 0, 1)
    )
  }

# MCMC settings
  ni <- 55000       # Number of draws from posterior (for each chain)
  nt <- 10          # Thinning rate
  nb <- 15000       # Number of draws to discard as burn-in
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
  out <- c(ests, mean(sw), mean(sbeta0), mean(sbetaT), mean(sk), mean(st0))
  names(out)[(length(out)-4):length(out)] <- c(
    'sw', 'sbeta0', 'sbetaT','sk','st0')  
  
# Return the result  
  return(list(
    out = out
    ))
  
} # End wrapper function
  
# Load libraries on workers -----
  sfLibrary(rjags)
  sfLibrary(rlecuyer)
  
# Start network random number generator -----
#  sfClusterSetupRNG()
  
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
  covres <- lapply(result, function(x) x[[c('out')]])
  covres <- data.frame(do.call(rbind, covres))
    
# Save results
  save(covres, file="covresult.rda")
  