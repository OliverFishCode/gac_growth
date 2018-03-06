# Front-end needs ---------------------------------------------------------
# Package install
  library(plyr)
  library(R2jags)
  # Need to install JAGS software before you can use this:
  # http://sourceforge.net/projects/mcmc-jags/files/

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
  nsamps = 30

# Create population of known ages with 30 individuals in each age class
  age = rep(seq(1, nages, 1), nsamps)
  
# Known parameters of VBGF for simulated population
  linf = 500
  k = 0.25
  t0 = -1.0
  sdlinf = 0.02
  sdk = 0.02
  sdt0 = 0.02
  
  slinf = rnorm(nages*nsamps, linf, sdlinf)    
  sk = rnorm(nages*nsamps, k, sdk)             
  st0 = rnorm(nages*nsamps, t0, sdt0) 
  sw = slinf*sk  
  
# Add a little random noise to each of the parameters
  slinf = linf + round(runif(nages*nsamps, -10, 10))
  sk = k + runif(nages*nsamps, -0.05, 0.05)
  st0 = t0 + runif(nages*nsamps, -1, 0.5)
  sw = slinf * sk
  
# Simulated length of individuals based on age and VBGF parameters
  slengq = (sw/sk)*(1-exp(-sk*(age-st0)))
  
# Put the data together in a dataframe  
  fish = data.frame(age, slinf, sk, st0, slengq)
  
# Model specification -----
# Write the model to a file
  modelString = "
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
    }"
  writeLines(modelString, con='vbModgq.txt')

# Model calibration -----
# Parameters monitored
  params = c('to', 'K', 'beta0')
  
# Package the data for JAGS
  vb_data = list(
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
  ni <- 55000      # Number of draws from posterior (for each chain)
  nt <- 10         # Thinning rate
  nb <- 15000      # Number of draws to discard as burn-in
  nc <- 3          # Number of chains

# Call jags and run the model
  vbModgq <- jags(data=vb_data, inits=inits, params, "vbModgq.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
    working.directory = getwd())

# Print a summary of the model
  print(vbModgq) 
  
# Model results -----
# Get posterior distributions for parameter estimates
  ek = vbModgq$BUGSoutput$sims.list$K
  et0 = vbModgq$BUGSoutput$sims.list$to
  ew = exp(vbModgq$BUGSoutput$sims.list$beta0)
  
# Check params for correlations
  summary(lm(ek~ew))
  plot(ek, exp(ew))
  
# Quick check of predictions
  ages=seq(1, nages, .5)
  Lt = (exp(mean(ew))/mean(ek))*(1-exp(-mean(ek)*(ages-mean(et0))))
  
  plot(age, slengq, pch=21, bg='black', cex=1.5)
  lines(ages, Lt, type='l', lty=1, lwd=2, col='blue')  
  
# Model predictions -----
# Growth curve - takes time to render
  # Predict Lt for each fish from model
    ages=seq(1, max(age), .5)
    preds = matrix(data = NA, nrow=length(ek), ncol=length(ages))
    for(i in 1:length(ek)){
      for(t in 1:length(ages)){
        preds[i, t] = (exp(ew[i])/ek[i])*
          (1-exp(-ek[i]*(ages[t]-et0[i])))
      }
    }
  # Make the posterior predictive plot  
    par(mar=c(4,4,1,1))
    plot(age, slengq, ylim=c(0, round(max(preds[ , ncol(preds)]), -1)),
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
    mtext(expression(paste(italic('w'))),
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

  # Omega
    hist(ew, col='gray87')                  # Estimate of omega
    abline(v=mean(sw), col = 'blue', lwd=2) # True mean
          
