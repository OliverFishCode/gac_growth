# Front-end needs ---------------------------------------------------------
# Package load
  library(plyr)
  library(R2jags)
  library(lubridate)

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
    
# Make a function for dividing that can be passed to apply
  divi <- function(x, y){
    out <- x/y
    return(out)
  }
  
# Make a function for subtraction that can be passed to apply
  subtract <- function(x, y){
    out <- x-y
    return(out)
  }
  
# Data read -----
# Read in the data
# These are data for sunfishes
# represented in the NYSDEC 2015
# State-wide fisheries database
  fish <- read.csv('sunnyAges.txt',
                   stringsAsFactors = FALSE)

# Remove everything but Pumpkinseed
  fish <- fish[fish$Name=='Pumpkinseed' & fish$Age <= 10,]
  
# Count the number of fish per waterbody
  tallies <-  data.frame(with(fish, table(Water, Name)))
  
# Remove tallies for waterbodies with fewer than 50 fish
  tallies <- tallies[tallies$Freq >= 50, ]
  #tallies <- tallies[tallies$Water!="Chautauqua Lake\\oiuy-098765",]
  
# Convert waterbody to character string for merging
  tallies$Water <- as.character(tallies$Water)
  
# Keep only those observations from the original data
# that have 50 or more samples
  fish <- fish[fish$Water %in% tallies$Water, ]
  
# Renumber row names
  row.names(fish) <- seq(1, nrow(fish), 1)
  
# Remove what appear to be erroneous age assignments or lengths
  # fish <- fish[-2183, ]
  # fish <- fish[-which(fish$Age==0 & fish$Length > 100),]
  # fish <- fish[-which(fish$Age==1 & fish$Length > 150),]
  # 
# Have a look at the length-at-age data
  plot(fish$Age, fish$Length)

# Calculate shoreline development index  
  fish$dev = (fish$ShrLen/.62)/(sqrt(4*pi*(fish$Sarea*.00156))) # Hutchinson (1957)
  
# Drop the two biggest lakes from the data set 
  fish <- fish[fish$Sarea < 10e3, ]
  
# Get year for each observation
  fish$Date <- as.Date(fish$Date, format="%m/%d/%Y")
  fish$year <- year(fish$Date)
  
# Get rid of missing values and age-0 fish
  fish <- fish[!is.na(fish$Age), ];
  fish <- fish[!is.na(fish$Length), ];
  fish <- fish[fish$Age!=0, ];  
  
# Create a new variable with combo of
# lake and region bc there are lakes
# in diff regions with same names.
  fish$WaterR <- paste0(fish$Water, fish$Region)
  
# Model specification -----
# . Write the model in bugs language ----
  modelString = "
    model{

        for(i in 1:N){
        # Likelihood
          Y[i] ~ dnorm(L[i], tau[Ti[i]])

        # Length described by Galluci & Quinn (1979)
          L[i] <- (w[i]/K[pop[i]])*(1-exp(-K[pop[i]]*(Ti[i]-to[pop[i]])))

        # Linear predictor of w
          log(w[i]) <- beta0[pop[i]] + betaY[region[i]]*dev[i]
        }

        # Priors
        for(j in 1:npops){
        # Priors on VBGF parameters (w defined below by linear model)
          # Brody growth coefficient
            lK[j] ~ dnorm(k.mu, k.tau)    
            logit(K[j]) <- lK[j]
      
          # Age at length zero
            to[j] ~ dnorm(t0.mu, t0.tau)T(-2,2)

        # Priors on parameters of linear model on w
          # Intercept
            beta0[j] ~ dnorm(b.mu, b.tau)
        }


        # Regional priors based on NYSDEC Regions
        for(r in 1:nregions){
          betaY[r] ~ dnorm(bd.mu, bd.tau)
        }
        

      # Global priors parameters above
        b.mu ~ dnorm(0, 0.001)
        b.tau ~ dgamma(0.01, 0.001)
        bd.mu ~ dnorm(0, 0.001)
        bd.tau ~ dgamma(0.01, 0.001)
        k.mu ~ dnorm(0, 0.001)
        k.tau ~ dgamma(0.01, 0.001)
        t0.mu ~ dnorm(0, 0.001)T(-2,0)
        t0.tau ~ dgamma(0.01, 0.001)

      # Prior distribution for precision at each age
      # This imposes a multiplicative error structure on length at age
        for(t in 1:Tmax){
          tau[t] ~ dgamma(0.01, 0.001)
        }
    }"

# Model calibration -----
# . Parameters monitored ----
  params = c('to', 'K', 'beta0', 'b.mu', 'b.tau', 'betaY',
             'bd.mu', 'bd.tau', 'k.mu', 'k.tau',
             't0.mu', 't0.tau')
  
# . Package the data for JAGS ----
  vb_data = list(
    Y = fish$Length,
    Ti = as.numeric(as.factor(fish$Age)),
    Tmax = max(as.numeric(as.factor(fish$Age))),
    N = nrow(fish),
    npops = length(unique(fish$WaterR)),
    pop = as.numeric(as.factor(fish$WaterR)),
    dev = as.vector(scale(fish$year)),
    nregions = length(unique(fish$Region)),
    region = as.numeric(as.factor(fish$Region)) 
  )
 
# . Initial values ---- 
  inits <- function(){
    list(
      b.mu = rnorm(1, 0, 1),
      b.tau = rgamma(1, 0.01, 1),
      bd.mu = rnorm(1, 0, 1),
      bd.tau = rgamma(1, 0.01, 1),
      k.mu = rnorm(1, 0, 1),
      k.tau = rgamma(1, 0.01, 1),
      t0.mu = runif(1, -2, 0),
      t0.tau = rgamma(1, 0.01, 1),      
      tau = rgamma(max(as.numeric(as.factor(fish$Age))), 0.01, 1)
    )
  }

# . MCMC settings ----
  ni <- 150000      # Number of draws from posterior (for each chain)
  nt <- 100         # Thinning rate
  nb <- 50000      # Number of draws to discard as burn-in
  nc <- 3          # Number of chains

# . Call jags and run the model ----
  vbModgq <- jags(data=vb_data, inits=inits, params,
                  textConnection(modelString),
                  n.chains = nc, n.thin = nt,
                  n.iter = ni, n.burnin = nb,
                  working.directory = getwd())

# Print a summary of the model
  print(vbModgq) 

# Save the results out to a file
  save(vbModgq, file='empiricalresult.rda')
  
# Model results -----
# . Get posteriors -----
# Read in the model file
  load('empiricalresult.rda')
  
# Get posterior distributions for parameter estimates
  ek <- vbModgq$BUGSoutput$sims.list$K
  et0 <- vbModgq$BUGSoutput$sims.list$to
  ew <- exp(vbModgq$BUGSoutput$sims.list$beta0)
  mu.w <- exp(vbModgq$BUGSoutput$sims.list$b.mu)
  betaY <- vbModgq$BUGSoutput$sims.list$betaY
  mu.k <- inv.logit(vbModgq$BUGSoutput$sims.list$k.mu)
  mu.t0 <- vbModgq$BUGSoutput$sims.list$t0.mu

# . State-wide growth curve -----
# Make a sequence of new ages for which we will predict lengths
  Age = seq(1, 13, 1)
  samp = nrow(ek)
  # Predict mean length at age for each sample
    preds = matrix(data = NA, nrow=samp, ncol=length(Age))
    for(i in 1:nrow(ek)){
      for(t in 1:length(Age)){
        preds[i, t] = mu.w[i]/mu.k[i]*(1-exp(-mu.k[i]*(Age[t]-mu.t0[i])))
      }
    }
  # Plot the raw data
    png(filename = 'statewide.png',
        height = 878, width = 1314,
        pointsize = 32)
    par(mar=c(5,5,1,1))
    plot(fish$Age[fish$Age!=0]+1,
         fish$Length[fish$Age!=0],
         ylim=c(0, 500),
         yaxt='n', xlab='Age (years)',
         ylab='Total length (mm)',
         xlim=c(0,12), axes = FALSE,
         pch = 21, bg=rgb(0.5,0.5,0.5,0.15),
         col=rgb(0.5,0.5,0.5,0.05),
         main='Statewide pumpkinseed growth')
  # Plot the posterior predictions
    for(i in 1:nrow(ek)){
      lines(x = Age, y = preds[i, ], col=rgb(.7,.7,.7,.05), lwd=1)
    }
  # Calculate the mean and 95% CRIs for posterior predictions
    muPred = apply(preds, 2, mean)
    lowPred = apply(preds, 2, low)
    upPred = apply(preds, 2, up)
  # Plot the mean and 95% CRI for predicted length at each age
    lines(Age, muPred, col='blue', lwd=2, lty=1)
    lines(Age, upPred, col='red', lwd=2, lty=2)
    lines(Age, lowPred, col='red', lwd=2, lty=2)
    axis(1, pos=0, at=seq(1,13,1), labels=seq(0,12,1))
    axis(2, pos=1, las=2)
    dev.off()  
    
# . Boxplot of k that shows 95% CRIs by waterbody ----
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ek,
                outline=FALSE, col='gray87', ylim=c(0,1),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ek,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,1),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(0,50,5), labels = seq(0,50,5))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic('k'['g'])), side=2, line=3.5)

# . Histogram of global omega (b.mu) ----
# Set plotting margins
  par(mar=c(5,5,1,1))
# Make the histogram
  hist(mu.w, main='', xlab=expression(mu[omega]),
       axes=FALSE, col = 'gray87', border='gray90',
       xlim = c(35,65), breaks=50)
# Add axes
  axis(side=1, pos=0)
  axis(side=2, pos=35, las=2)
# Add line for observed length at age 0
  abline(v=mean(mu.w), col='gray40', lwd=2, lty=3)
  
      
# . Boxplot of omega that shows 95% CRIs by waterbody ----
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ew,
                outline=FALSE, col='gray87', ylim=c(0,125),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ew,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,125),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(0,50,5), labels = seq(0,50,5))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(omega['g']), side=2, line=3.5)
  abline(h=median(ew))
  
# . Boxplot of coefficient that shows 95% CRIs by waterbody ----
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(betaY,
                outline=FALSE, col='gray87', ylim=c(-.5,.5),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    betaY,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(-.5,.5),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', boxwex=0.3)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=c(1,seq(2,8,1)), labels = c(1,seq(3,9,1)))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='NYSDEC management region', side=1, line=3.5)
  mtext(text=expression(beta['g']), side=2, line=3.5)
  abline(h=0, lty=2)
  

# . Compare omegas to state-wide average ----
  # Calculate ratio of waterbody omega to statewide average omega
  reg <-data.frame(unique(cbind(fish$WaterR, fish$Region)))
  regs <- reg[with(reg, order(X1)),]
  row.names(regs) <- seq(1:nrow(regs))
  regs2 <- regs[with(regs, order(X2)),]
  reg <- as.factor(regs2[,2])
  
  quotient <- apply(ew, 2, divi, mu.w)
  quotient <- quotient[,as.numeric(row.names(regs))]
  #diff <- apply(ew, 2, subtract, mu.w)
  
  # Set up plotting file
    png(filename = 'deltas.png',
      height = 878, width = 1314,
      pointsize = 32)
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(quotient,
                outline=FALSE, col='gray87', ylim=c(.5,2),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    quotient,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(.5,2), xlim=c(0,51),
      boxfill=gray.colors(length(unique(reg)))[reg],
      col.axis='white',
      staplewex=0, 
      staplecol='gray40',
      whisklty=1,
      whiskcol='gray40',
      whisklwd=2, 
      boxcol='gray40',
      boxlwd=1,
      medcol='gray40',
      xaxt='n')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(0,50,5), labels = seq(0,50,5))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match

  # Add axis labels
  mtext(text='Stock ID', side=1, line=3.5)
  mtext(text=expression(delta['g']), side=2, line=3.5)
  #abline(h=1, lty=2)  
  abline(v= which(!duplicated(reg))[-1]-.5,
         col='gray30')
  for(i in 1:length(unique(reg))){
    text(x=c((which(!duplicated(reg))[-1]-.5),52)[i],
         y=2,
         labels=unique(reg)[i],
         adj=c(rep(1.2, 3),1,rep(1.2, 3),.5)[i])
  }
  dev.off()
# . Write lake-specific quotients to a file -----
  quoti <- apply(quotient, 2, mean)
  outquote <- data.frame(regs, quoti)
  boxplot(quoti~X2, data=outquote)
  
# . Figure 4 ----
# .. Load empirical results ----
  load('empiricalresult.rda')
  
# Get posterior distributions for parameter estimates
  ek <- vbModgq$BUGSoutput$sims.list$K
  et0 <- vbModgq$BUGSoutput$sims.list$to
  ew <- exp(vbModgq$BUGSoutput$sims.list$beta0)
  mu.w <- exp(vbModgq$BUGSoutput$sims.list$b.mu)
  betaY <- vbModgq$BUGSoutput$sims.list$betaY
  mu.k <- inv.logit(vbModgq$BUGSoutput$sims.list$k.mu)
  mu.t0 <- vbModgq$BUGSoutput$sims.list$t0.mu  

# .. Set up image file ----
  tiff(filename = 'Figure4.tif',
      height = 2000, width = 2600,
      res=350, pointsize = 10)  
  
# Set up graphical parameters
  par(mfrow=c(2,2), oma=c(1,1,0,0))
  
# .. State-wide growth curve ----- 
# Make a sequence of new ages for which we will predict lengths
  Age = seq(1, 13, 1)
  samp = nrow(ek)
  # Predict mean length at age for each sample
    preds = matrix(data = NA, nrow=samp, ncol=length(Age))
    for(i in 1:nrow(ek)){
      for(t in 1:length(Age)){
        preds[i, t] = mu.w[i]/mu.k[i]*(1-exp(-mu.k[i]*(Age[t]-mu.t0[i])))
      }
    }
  # Plot the raw data
    par(mar=c(4.5,1.25,1.5,2))
    plot(fish$Age[fish$Age!=0]+1,
         fish$Length[fish$Age!=0],
         ylim=c(0, 300),
         yaxt='n', xlab='',
         ylab='', xaxt='n',
         xlim=c(0,12), axes = FALSE,
         pch = 21, bg=rgb(0.5,0.5,0.5,0.15),
         col=rgb(0.5,0.5,0.5,0.05),
         main='')
  # Plot the posterior predictions
    for(i in 1:nrow(ek)){
      lines(x = Age, y = preds[i, ], col=rgb(.7,.7,.7,.05), lwd=1)
    }
  # Calculate the mean and 95% CRIs for posterior predictions
    muPred = apply(preds, 2, mean)
    lowPred = apply(preds, 2, low)
    upPred = apply(preds, 2, up)
  # Plot the mean and 95% CRI for predicted length at each age
    lines(Age, muPred, col='black', lwd=1, lty=1)
    lines(Age, upPred, col='black', lwd=1, lty=2)
    lines(Age, lowPred, col='black', lwd=1, lty=2)
  # Add axes and labels
    axis(1, pos=0, at=seq(1,13,1), labels=seq(0,12,1))
    axis(2, pos=1, las=2)
    mtext('Age (years)', side=1, line=3)
    par(xpd=NA)
    text(x=-1.25, y=150, 'Total length (mm)', srt=90)
    #mtext('Total length (mm)', outer=TRUE, line=0)

    
# .. Boxplot of coefficient that shows 95% CRIs by waterbody ----
  # Set margins
  par(mar=c(5,4,2,2), xpd=FALSE)
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(betaY,
                outline=FALSE, col='gray87', ylim=c(-.5,.5),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    betaY,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(-.5,.5),
      boxfill='gray87', col.axis='white',
      staplewex=0, staplecol='gray40', whisklty=1, whiskcol='gray40',
      whisklwd=1, boxcol='gray40', boxlwd=1,
      medcol='gray40', boxwex=0.3)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=c(1,seq(2,8,1)), labels = c(1,seq(3,9,1)))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Management region', side=1, line=3.5)
  mtext(text=expression(beta['g']), side=2, line=3.5)
  abline(h=0, lty=2)
    
# .. Boxplot of omega that shows 95% CRIs by waterbody ----
  # Calculate ratio of waterbody omega to statewide average omega
  reg <-data.frame(unique(cbind(fish$WaterR, fish$Region)))
  regs <- reg[with(reg, order(X1)),]
  row.names(regs) <- seq(1:nrow(regs))
  regs2 <- regs[with(regs, order(X2)),]
  reg <- as.factor(regs2[,2])

  ew <- ew[,as.numeric(row.names(regs))]  
  
  # Set margins
  par(mar=c(5,4,2,2))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ew,
                outline=FALSE, col='gray87', ylim=c(0,125),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ew,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,125),
      boxfill=gray.colors(length(unique(reg)))[reg],
      col.axis='white',
      staplewex=0, staplecol='gray40', whisklty=1, whiskcol='gray40',
      whisklwd=1, boxcol='gray40', boxlwd=1,
      medcol='gray40', xaxt='n')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(0,50,5), labels = seq(0,50,5))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Stock ID', side=1, line=3.5)
  mtext(text=expression(omega['g']), side=2, line=3.5)
  #abline(h=median(ew))
    abline(v= which(!duplicated(reg))[-1]-.5,
         col='gray30')
  for(i in 1:length(unique(reg))){
    text(x=c((which(!duplicated(reg))[-1]-.5),52)[i],
         y=125,
         labels=unique(reg)[i],
         adj=c(rep(1.2, 3),1,rep(1.2, 3),.5)[i])
  }

# .. Compare omegas to state-wide average ----
  # Calculate ratio of waterbody omega to statewide average omega
  reg <-data.frame(unique(cbind(fish$WaterR, fish$Region)))
  regs <- reg[with(reg, order(X1)),]
  row.names(regs) <- seq(1:nrow(regs))
  regs2 <- regs[with(regs, order(X2)),]
  reg <- as.factor(regs2[,2])
  
  quotient <- apply(ew, 2, divi, mu.w)
  quotient <- quotient[,as.numeric(row.names(regs))]
  
  # Set margins
  par(mar=c(5,4,2,2))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(quotient,
                outline=FALSE, col='gray87', ylim=c(.5,2),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    quotient,
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(.5,2), xlim=c(0,51),
      boxfill=gray.colors(length(unique(reg)))[reg],
      col.axis='white',
      staplewex=0, 
      staplecol='gray40',
      whisklty=1,
      whiskcol='gray40',
      whisklwd=1, 
      boxcol='gray40',
      boxlwd=1,
      medcol='gray40',
      xaxt='n')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(0,50,5), labels = seq(0,50,5))
  axis(side=2, las=2)
  # Add axis labels
  mtext(text='Stock ID', side=1, line=3.5)
  mtext(text=expression(paste(Delta, omega['g'])), side=2, line=3.5)
  abline(v= which(!duplicated(reg))[-1]-.5,
         col='gray30')
  for(i in 1:length(unique(reg))){
    text(x=c((which(!duplicated(reg))[-1]-.5),52)[i],
         y=2,
         labels=unique(reg)[i],
         adj=c(rep(1.2, 3),1,rep(1.2, 3),.5)[i])
  }    
    
# .. Close graphics device ----
  dev.off()  
    
  
  