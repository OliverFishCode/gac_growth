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
  linf = vector(mode='list', length=npops)
  k = vector(mode='list', length=npops)
  t0 = vector(mode='list', length=npops)

  for(i in 1:npops){
    linf[[i]] = rep(seq(550, 350, -200/npops)[i], nsamps*nages)
    k[[i]] = rep(seq(0.1, 0.5, 0.40/npops)[i], nsamps*nages)
    t0[[i]] = rep(seq(-3, 1, -2/npops)[i], nsamps*nages)

  }
  sdlinf = 20
  sdk = 0.02
  sdt0 = 0.02  
  
  pops = sort(rep(seq(1,10, 1), nsamps*nages))
  ages = rep(seq(1,10, 1), nsamps*nages)
  
  parms = data.frame(pops, ages,
                     unlist(linf), 
                     unlist(k),
                     t0, sdlinf, sdk, sdt0)
  names(parms) = c('pops', 'ages', 'linf', 'k', 't0', 'sdlinf', 'sdk', 'sdt0')
  
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
  ni <- 5500     # Number of draws from posterior (for each chain)
  nt <- 10        # Thinning rate
  nb <- 1500     # Number of draws to discard as burn-in
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
  ew = vbModgq$BUGSoutput$sims.list$beta0
  
# # Check params for correlations
#   summary(lm(ek~ew))
#   plot(ek, exp(ew))
#   
# # Quick check of predictions
#   ages=seq(1, nages, .5)
#   Lt = (exp(mean(ew))/mean(ek))*(1-exp(-mean(ek)*(ages-mean(et0))))
#   
#   plot(age, slengq, pch=21, bg='black', cex=1.5)
#   lines(ages, Lt, type='l', lty=1, lwd=2, col='blue')  
#   
# # Model predictions -----
# # Growth curve - takes time to render
#   # Predict Lt for each fish from model
#     ages=seq(1, max(age), .5)
#     preds = matrix(data = NA, nrow=length(ek), ncol=length(ages))
#     for(i in 1:length(ek)){
#       for(t in 1:length(ages)){
#         preds[i, t] = (exp(ew[i])/ek[i])*
#           (1-exp(-ek[i]*(ages[t]-et0[i])))
#       }
#     }
#   # Make the posterior predictive plot  
#     par(mar=c(4,4,1,1))
#     plot(age, slengq, ylim=c(0, round(max(preds[ , ncol(preds)]), -1)),
#          yaxt='n',
#          xlab='',
#          ylab='',
#          xlim=c(0, max(ages)),
#          axes = FALSE,
#          pch = 21, bg='black', col='black', main='')
#   # Plot the posterior predictions
#     for(i in 1:length(ek)){
#       lines(x = ages, y = preds[i, ], col=rgb(.7,.7,.7,.05), lwd=1)
#     }
#   # Calculate the mean and 95% CRIs for posterior predictions
#     muPred = apply(preds, 2, mean)
#     lowPred = apply(preds, 2, low)
#     upPred = apply(preds, 2, up)
#   # Plot the mean and 95% CRI for predicted length at each age
#     lines(ages, muPred, col='blue', lwd=2, lty=1)
#     lines(ages, upPred, col='red', lwd=2, lty=2)
#     lines(ages, lowPred, col='red', lwd=2, lty=2)
#     axis(1, pos=0)
#     axis(2, pos=0, las=2)  
#     mtext(expression(paste('Age (years)')),
#           side=1, line=2.5)
#     mtext(expression(paste(italic('w'))),
#           side=2, line=2.5)  
#   