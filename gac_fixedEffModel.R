# Front-end needs ---------------------------------------------------------
# Package load
  library(plyr)
  library(R2jags)

# Make function for inverting logit
  inv.logit=function(x){
    exp(x)/(1+exp(x))
  }

# Make a function to get lower 95% credible limit with short name
  low = function(x){
    quantile(x, probs=c(0.025))
  }

# Make a function to get upper 95% credible limit with short name
  up = function(x){
    quantile(x, probs=c(0.975))
  }
    
# Data definition -----
# Number of age classes
  nages = 10
  
# Number of samples in each age class
  nsamps = 10

# Create population of known ages with 30 individuals in each age class
  age = rep(seq(1, nages, 1), nsamps)
  
# Create a known number of populations over which to replicate the 
# age structures
  npops = 10
  
# Known parameters of VBGF for simulated populations, with linf and
# k being pop-specific
  # Pre-allocate vectors to hold known von Bertalanffy params
    linf = vector(mode='list', length=npops)
    k = vector(mode='list', length=npops)
    t0 = vector(mode='list', length=npops)
  
  # Assign known values of von Bertalanffy params
  # using a sequence of equally spaced values for each
    for(i in 1:npops){
      linf[[i]] = rep(seq(550, 350, -200/npops)[i], nsamps*nages)
      k[[i]] = rep(seq(0.1, 0.5, 0.40/npops)[i], nsamps*nages)
      t0[[i]] = rep(seq(-3, 1, 2/npops)[i], nsamps*nages)
    }
  
  # Define standard deviations for each parameter that
  # will be used to draw simulated parameter values
    sdlinf = 20
    sdk = 0.02
    sdt0 = 0.02  
  
  # Make sequences for population number and age for
  # each of the simulated fish
    pops = sort(rep(seq(1,10, 1), nsamps*nages))
    ages = rep(seq(1,10, 1), nsamps*length(unique(pops)))
  
  # Put the known parameters in a dataframe  
    parms = data.frame(pops, ages,
                       unlist(linf), 
                       unlist(k),
                       unlist(t0),
                       sdlinf,
                       sdk,
                       sdt0)
  # Give it some manageable names  
    names(parms) = c('pops',
                     'ages',
                     'linf',
                     'k',
                     't0',
                     'sdlinf',
                     'sdk',
                     'sdt0')
    
  # Simulate nsamps number of von Bert parameters for each age in each
  # population
    slinf =  rnorm(nrow(parms), parms$linf, parms$sdlinf)             
    sk =  rnorm(nrow(parms), parms$k, parms$sdk)             
    st0 = rnorm(nrow(parms), parms$t0, parms$sdt0)       
    sw = slinf * sk 
    
  # Simulated length of individuals based on age and VBGF parameters
    slengq = (sw/sk)*(1-exp(-sk*(age-st0)))
    
  # Put the data together in a dataframe  
    fish = data.frame(parms, slinf, sk, st0, slengq)
  
# Model specification -----
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
            beta0[j] ~ dnorm(0, 0.001)
        }

      # Prior distribution for precision at each age
      # This imposes a multiplicative error structure on length at age
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.01, 0.001)
        }
    }"
  writeLines(modelString, con='vbModgq_raneff.txt')

# Model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0')
  
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
      beta0 = rnorm(length(unique(fish$pops)), 0, 1),     
      K = runif(length(unique(fish$pops)), 0, 1),
      to = runif(length(unique(fish$pops)), -10, 1),
      tau = rgamma(max(fish$age), .01, 1)
    )
  }

# MCMC settings
  ni <- 5500      # Number of draws from posterior (for each chain)
  nt <- 10        # Thinning rate
  nb <- 1500      # Number of draws to discard as burn-in
  nc <- 3         # Number of chains

# Call jags and run the model
  vbModgq <- jags(data=vb_data, inits=inits, params, "vbModgq_raneff.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
    working.directory = getwd())

# Print a summary of the model
  print(vbModgq) 

# Model results -----
# Get posterior distributions for parameter estimates
  ek = vbModgq$BUGSoutput$sims.list$K
  et0 = vbModgq$BUGSoutput$sims.list$to
  ew = exp(vbModgq$BUGSoutput$sims.list$beta0)
    
# Make some quick boxplots to make sure the
# posteriors (boxes) follow simulated parameter values (points)
  # Graphics window margins
    par(mar=c(5,5,1,1))
  # k
    # Posterior distribution for each pop    
      boxplot(ek, col='gray87', outline=FALSE,
              xlab='Population (i)',
              ylab=expression(italic('k')[italic('i')]),
              axes=FALSE,
              ylim=c(0,1), xlim=c(0,11)
              )  
    # Axes  
      axis(1, pos=0)
      axis(2, las=2, pos=0)
    # Known means for each pop
      points(x=unique(pops), 
             y=unique(unlist(k)),
             pch=21, cex=1.5, bg='blue', col='blue')
      
  # Omega
    # Posterior distribution for each pop
      boxplot(ew, col='gray87', outline=FALSE,
              xlab='Population (i)', ylab=expression(omega[italic('i')]),
              axes=FALSE,
              ylim=c(0,300), xlim=c(0,11)
              )  
    # Axes
      axis(1, pos=0)
      axis(2, las=2, pos=0)
    # Known means for each pop
      points(x=unique(pops), 
         y=unique(unlist(k))*unique(unlist(linf)),
         pch=21, cex=1.5, bg='red', col='red')
    