# Front-end needs ---------------------------------------------------------
# Package install
  library(plyr)
  library(R2jags)
  # Need to install JAGS software before you can use this:
  # http://sourceforge.net/projects/mcmc-jags/files/

# Set seed for RNGesus
  set.seed(607)

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
  nsamps = 20

# Create population of known ages with 30 individuals in each age class
  age = rep(seq(1, nages, 1), nsamps)
  
# Known parameters of VBGF for simulated population
  k = 0.35      # Mean of k
  sdk = 0.02    # SD of k
  t0 = -1.0     # Mean of t0
  sdt0 = 0.02   # SD of t0
  beta0 = 5     # Mean of intercept for glm on omega
  sdbeta0 =.2   # SD of intercept
  betaT = .12   # Mean effect of temperature
  sdbetaT = .02 # SD of temperature effect
  temp = runif(length(age), 15, 25) # Simulated temperatures
  
# Add some error to each of the parameters by 
# drawing each from distributions
  sk = rnorm(nages*nsamps, k, sdk)             
  st0 = rnorm(nages*nsamps, t0, sdt0)           
  sbeta0 = rnorm(nages*nsamps, beta0, sdbeta0) 
  sbetaT = rnorm(nages*nsamps, betaT, sdbetaT)
  sw = exp(sbeta0 + sbetaT*as.vector(scale(temp)))   
  
# Simulated length of individuals based on age and VBGF parameters
  slengq = (sw/sk)*(1-exp(-sk*(age-st0))) # Galluci and Quinn (1979)

# Put the simulated data together in a dataframe
  fish = data.frame(age,    # Fish age
                    temp,   # Temperature
                    sw,     # L-infinity
                    sk,     # Brody growth coefficient
                    st0,    # t0
                    slengq  # Length from gq params and age- identical
                    )
  
# Covariate model specification -----
# Write the model to a file
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
  writeLines(modelString, con='vbModgq.txt')

# Covariate model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0', 'betaT', 'tau')
  
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
  ni = 15000     # Number of draws from posterior (for each chain)
  nt = 5         # Thinning rate
  nb = 5000      # Number of draws to discard as burn-in
  nc = 3         # Number of chains

# Call jags and run the model
  vbModgq <- jags(data=vb_data, inits=inits, params, "vbModgq.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
    working.directory = getwd())

# Covariate model results -----
# Print a summary of the model
  print(vbModgq) 
  
# Get posterior distributions for parameter estimates
# Posterior samples from each MCMC run are elements of 
# sims.list. Can be vectors or matrices depending on 
# model formulation. Here, the prefix 'e' stands for 'estimate'
  ek = vbModgq$BUGSoutput$sims.list$K
  et0 = vbModgq$BUGSoutput$sims.list$to
  ebeta0 = vbModgq$BUGSoutput$sims.list$beta0
  ebetaT = vbModgq$BUGSoutput$sims.list$betaT
  ew = exp(ebeta0 + ebetaT*mean(as.vector(scale(fish$temp))))

# Check params for correlations.
  summary(lm(ek~ew))
  plot(ek, ew)
  
# Quick check of predictions
  # Predict mean length at age using coefficient estimates
    ages = seq(1, nages, 1)
    Lt = (mean(ew)/mean(ek))*(1-exp(-mean(ek)*(ages-mean(et0))))
  # Plot the raw data followed by a blue line for mean prediction
    plot(age, slengq, pch=21, bg='black', cex=1.5)
    lines(ages, Lt, type='l', lty=1, lwd=2, col='blue')  
  
# Model predictions -----
    
# **WARNING** These plots take a while to render
    
# Growth curve
  # Predict Lt for each fish from model coefficients
    ages=seq(1, max(age), .5)
    preds = matrix(data = NA, nrow=length(ek), ncol=length(ages))
    for(i in 1:length(ek)){
      for(t in 1:length(ages)){
        preds[i, t] = (ew[i]/ek[i])*
          (1-exp(-ek[i]*(ages[t]-et0[i])))
      }
    }
  # Make the posterior predictive plot  
    par(mar=c(4,4,1,1))
    plot(age, slengq, ylim=c(0, round(max(fish$slengq))),
         yaxt='n',
         xlab='',
         ylab='',
         xlim=c(0, max(ages)),
         axes = FALSE,
         pch = 21, bg='black', col='black', main='')
  # Plot the posterior predictions
    for(i in 1:length(ek)){
      lines(x = ages, y = preds[i, ], col=rgb(.7,.7,.7,.05), lwd=1)
    }
  # Calculate the mean and 95% CRIs for posterior predictions
    muPred = apply(preds, 2, mean)
    lowPred = apply(preds, 2, low)
    upPred = apply(preds, 2, up)
  # Plot the mean and 95% CRI for predicted length at each age
    lines(ages, muPred, col='blue', lwd=2, lty=1)
    lines(ages, upPred, col='red', lwd=2, lty=2)
    lines(ages, lowPred, col='red', lwd=2, lty=2)
    axis(1, pos=0)
    axis(2, pos=0, las=2)  
    mtext(expression(paste('Age (years)')),
          side=1, line=2.5)
    mtext(expression(paste('Length (mm)')),
          side=2, line=2.5)  
    
# Checks for accuracy -----
# Parameter recovery comparisons. How do posteriors match up with
# known values used for simulation? Plot the posterior 
# distribution for each parameter of interest against the mean
# used for simulation.
  # k
    hist(ek, col='gray87')
    abline(v=mean(sk), col = 'blue', lwd=2) # True mean

  # t0
    hist(et0, col='gray87')
    abline(v=mean(st0), col = 'blue', lwd=2) # True mean

  # Effect of temperature
    hist(ebetaT, col='gray87')
    abline(v=mean(sbetaT), col = 'blue', lwd=2) # True mean

  # Omega
    hist(ew, col='gray87')                  # Estimate of omega
    abline(v=mean(sw), col = 'blue', lwd=2) # True mean
    
    
    
    
  
