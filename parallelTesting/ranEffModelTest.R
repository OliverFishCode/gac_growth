# Front-end needs ---------------------------------------------------------
# Load R packages
  library(rjags)
  library(snowfall)
  library(rlecuyer)

# Get number of cores
  args = commandArgs(trailingOnly = TRUE);
  ncpus = args[1];
  #ncpus = 3

# Initialize snowfall
  sfInit(parallel = TRUE, cpus=ncpus, type="SOCK")

# Wrapper fxn -----
# Define the wrapper function to call in parallel,
# which is really just the whole simulation
  wrapper <- function(x){
    
# . Data definition -----
# Number of age classes
  nages <- 10
  
# Number of samples in each age class
  nsamps <- 30

# Create population of known ages with 30 individuals in each age class
  age <- rep(seq(1, nages, 1), nsamps)
  
# Create a known number of populations over which to replicate the 
# age structures
  npops <- 10
  
# Known parameters of VBGF for simulated populations, with linf and
# k being pop-specific
  # Pre-allocate vectors to hold known von Bertalanffy params
    sw <- vector(mode='list', length=npops)
    sk <- vector(mode='list', length=npops)
    st0 <- vector(mode='list', length=npops)
  
  # Assign known values of von Bertalanffy params
  # using a sequence of equally spaced values for each
    for(i in 1:npops){
      sw[[i]] <- rep(seq(50, 250, 200/npops)[i], nsamps*nages)
      sk[[i]] <- rep(seq(0.1, 0.5, 0.40/npops)[i], nsamps*nages)
      st0[[i]] <- rep(-1, nsamps*nages)
    }
    sw <- unlist(sw)
    sk <- unlist(sk)
    st0 <- unlist(st0)
    

  # Make sequences for population number and age for
  # each of the simulated fish
    pops <- sort(rep(seq(1,10, 1), nsamps*nages))
    ages <- rep(seq(1,10, 1), nsamps*length(unique(pops)))
  
  # Put the known parameters in a dataframe
    parms <- data.frame(pops, ages,
                       unlist(sw),
                       unlist(sk),
                       unlist(st0))
  # # Give it some manageable names  
    names(parms) <- c('pops',
                     'ages',
                     'sw',
                     'sk',
                     'st0')
    
  # Simulate length based on age and vonB params
  slengq <- 
  rnorm(n=nages*nsamps*npops,
        mean=(sw/sk)*(1-exp(-sk*(ages-st0))),
        sd=10)
  
  # Put the data together in a dataframe  
    fish <- data.frame(parms, slengq)
  
# . Model specification -----
# Write the model to a file
  modelString = "
    model{

        for(i in 1:N){
        # Likelihood
          Y[i] ~ dnorm(L[i], tau[Ti[i]])

        # Length described by Galluci & Quinn (1979)
          L[i] <- (w[i]/K[pop[i]])*(1-exp(-K[pop[i]]*(Ti[i]-to[pop[i]])))

        # Linear predictor of w
          log(w[i]) <- beta0[pop[i]]
        }

        for(j in 1:npops){
        # Priors on VBGF parameters (w defined below by linear model)
          # Brody growth coefficient
            K[j] ~ dunif(0, 1)
          # Age at length zero
            to[j] ~ dunif(-10, 1)

        # Priors on parameters of linear model on w
          # Intercept
            beta0[j] ~ dnorm(b.mu, b.tau)
        }

      # Global priors on mu and tau for beta0[j] above
        b.mu ~ dnorm(0, 0.001)
        b.tau ~ dgamma(0.01, 0.001)

      # Prior distribution for precision at each age
      # This imposes a multiplicative error structure on length at age
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.01, 0.001)
        }
    }"

# Model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0', 'b.mu', 'b.tau')
  
# Package the data for JAGS
  vb_data = list(
    Y = fish$slengq,
    Ti = fish$ages,
    Tmax = max(fish$ages),
    N = nrow(fish),
    npops = length(unique(fish$pops)),
    pop = fish$pops
  )
 
# Initial values 
  inits <- function(){
    list(
      b.mu = rnorm(1, 0, 1),
      b.tau = rgamma(1, 0.01, 1),
      K = runif(length(unique(fish$pops)), 0, 1),
      to = runif(length(unique(fish$pops)), -10, 1),
      tau = rgamma(max(fish$age), 0.01, 1)
    )
  }

# MCMC settings
  ni <- 5500      # Number of draws from posterior (for each chain)
  nt <- 10        # Thinning rate
  nb <- 1500      # Number of draws to discard as burn-in
  nc <- 3         # Number of chains

# Call jags and run the model, re-run if crashes due to
# bad initial values
  fin <- NULL
  attempt <- 0
  while(is.null(fin)){
    attempt <- attempt + 1
    try({
      vbMod <- jags.model(file=textConnection(modelString) , data=vb_data,
                            inits=inits, n.adapt = 100, quiet=TRUE)
      update(vbMod, nb)
      fin <- coda.samples(vbMod, params, ni)
      }

    )
  }

# Save the estimates and the known parameter values
# Save the estimates and the known parameter values
  ests <- summary(fin)[[1]][,1]
  out <- c(ests,
           unique(sw),#unlist(lapply(sw, unique)),
           unique(sk),#unlist(lapply(sk, unique)),
           rep(unique(st0), 10)#unlist(lapply(st0, unique))
           )
  names(out)[(length(out)-29):length(out)] <- c(
    paste('sw', seq(1,10,1), sep=''),
    paste('sk', seq(1,10,1), sep=''),
    paste('st0', seq(1,10,1), sep=''))  
  
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
  ranres <- lapply(result, function(x) x[[c('out')]])
  ranres <- data.frame(do.call(rbind, ranres))
    
# Save results
  save(ranres, file="ranresult.rda")
  