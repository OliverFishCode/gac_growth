# Front-end needs ----
# Function definitions
# Make a function to get lower 95% credible limit
# with short name
  low = function(x){
    quantile(x, probs=c(0.025))
  }

# Make a function to get upper 95% credible limit
# with short name
  up = function(x){
    quantile(x, probs=c(0.975))
  }

# Base model -----
load('result.rda')
head(res)
nrow(res)

# . Fig. 1 -----
  # Set up graphics device
  tiff(
    filename = "Figure1.tif",
    width = 1460,
    height = 970,
    pointsize = 8,
    res = 350
  )

  # Estimation accuracy for k
    par(mfrow=c(1, 3), oma=c(1,1,1,1), mar=c(4,4,1,1))
    hist(res$K, col='gray87', xlab='K',
         ylab = 'Frequency', yaxt='n', xaxt='n',
         xlim = c(0.2, .4), main='',
         border='gray87', ylim=c(0,250),
         breaks = 15, cex.lab=1.25
         )
    abline(v=mean(res$sk), col='black', lty=2, lwd=1)
    axis(side=1, pos=0)    
    axis(side=2, pos=.2, las=TRUE)
  
  # Estimation accuracy for omega
    par(mar=c(4,2,1,1))
    hist(exp(res$beta0), col='gray87', 
         xlab=expression(omega),
         ylab = 'Frequency', yaxt='n', xaxt='n',
         xlim = c(50, 150), main='',
         border='gray87', ylim=c(0,250),
         breaks = 15, cex.lab=1.25
         )
    abline(v=mean(res$sw), col='black',
           lty=2, lwd=1)
    axis(side=1, pos=0)    

  # Estimation accuracy for t0
    par(mar=c(4,2,1,1))
    hist(res$to, col='gray87', 
         xlab=expression('t'['0']),
         ylab = 'Frequency', yaxt='n', xaxt='n',
         xlim = c(-2, 0), main='',
         border='gray87', ylim=c(0,250),
         breaks = 30, cex.lab=1.25
         )
    abline(v=mean(res$st0), col='black', lty=2, lwd=1)
    axis(side=1, pos=0)    
    dev.off()

  # Accuracy and precision
    # k
    mean(res$K)
    quantile(res$K, probs = c(0.025, 0.975))     
    # omega
    mean(exp(res$beta0))
    quantile(exp(res$beta0), probs = c(0.025, 0.975)) 
    # t0
    mean(res$to)
    quantile(res$to, probs = c(0.025, 0.975))      
    
# Covariate model -----
  load('covresult.rda')
  head(covres)
  nrow(covres)

# Estimation accuracy for k
  par(mar=c(5,5,1,1))
  hist(covres$K, col='gray87', xlab=expression(italic('k')),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(.2, .4), main='')
  abline(v=mean(covres$sk), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=.2, las=TRUE)

# Estimation accuracy for betaT
  par(mar=c(5,5,1,1))
  hist(covres$betaT, col='gray87', 
       xlab=expression(beta['T']),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(.1, .14), main='')
  abline(v=mean(covres$sbetaT), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=.1, las=TRUE)

# Estimation accuracy for w
# Calculate w as function of linear predictor
  ew = exp(covres$beta0)
# Make plot
  par(mar=c(5,5,1,1))
  hist(ew, col='gray87', 
       xlab=expression(omega),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(0, 300), main='')
  abline(v=mean(covres$sw), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=0, las=TRUE)

# Data visualization for covariate effect on omega
  # Make new sequence of standardized temperatures
  newT <- seq(-2, 2, .01)
  # Make empty matrix to hold predictions
  pred <- matrix(NA, nrow=nrow(covres), ncol=length(newT))
  # Make mean prediction for each model at each temp
  for(i in 1:nrow(pred)){
    for(t in 1:ncol(pred)){
      pred[i,t] <- exp(covres$beta0[i] + covres$betaT[i]*newT[t])
    }
  }
  
  # Make the plot
  # Set graphical params
  par(mar=c(5,5,1,1))
  # Plot the first row of prediction against covariate
  plot(newT, pred[1, ], type='l', lty=1, lwd=1, 
       col=rgb(.4,.4,.4,.1), yaxt='n', xaxt='n',
       xlab='Covariate', ylab=expression(omega),
       ylim=c(110, 210))
  # Add remaining predictions
  for(i in 1:nrow(pred)){
    lines(newT, pred[i, ], lty=1, lwd=1, 
       col=rgb(.4,.4,.4,.1) )
  }
  # Add means and CIs
  lines(newT, apply(pred, 2, mean), lty=1, lwd=2, col='black')
  lines(newT, apply(pred, 2, up), lty=2, lwd=1, col='black')
  lines(newT, apply(pred, 2, low), lty=2, lwd=1, col='black')
  # Add axes
  axis(side=1)    
  axis(side=2,las=2)    
  
# Estimation accuracy for t0
# Make plot
  par(mar=c(5,5,1,1))
  hist(covres$to, col='gray87', 
       xlab=expression(t[0]),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(-2, 0), main='')
  abline(v=mean(covres$st0), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=-2, las=TRUE)  
  
  
# Fixed- and random-effects models -----
  load('fixedresult.rda')
  head(fixedres)
  nrow(fixedres)
  
  load('ranresult.rda')
  head(ranres)
  nrow(ranres)  

# . Fig. 3 ----
  # Set up graphics device
  tiff(
    filename = "Figure3.tif",
    width = 1280,
    height = 1280,
    pointsize = 8,
    res = 350
  )

 # Set up plotting window
    par(mfrow=c(3, 2), oma=c(4,5,2,1), mar=c(1,1,1,1))
 # Estimation accuracy for k from fixed eff model
 # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(fixedres[ , grep(pattern="K.", x=names(fixedres))],
                outline=FALSE, col='gray87', ylim=c(0,.5),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    fixedres[ , grep(pattern="K.", x=names(fixedres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,.5),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = FALSE)
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text=expression('K'[g]), side=2, line=3.5)
  # Known parameter values
  points(1:10,
         fixedres[1 , grep(pattern="sk", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )
  # Panel label
  mtext(text='Fixed effects model', side=3, line=1)

  
 # Estimation accuracy for k from random effects model
 # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ranres[ , grep(pattern="K.", x=names(ranres))],
                outline=FALSE, col='gray87', ylim=c(0,.5),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ranres[ , grep(pattern="K.", x=names(ranres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,.5),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = FALSE)
  axis(side=2, las=2, labels=FALSE)
  # Close it up with a box so axis styles match
  box()
  # Known values
  points(1:10,
         ranres[1 , grep(pattern="sk", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )
  # Panel label
  mtext(text='Random effects model', side=3, line=1)

  
  
# Estimation accuracy for omega from fixed effects model
# Compare estimated omega to true values using 
# a boxplot that shows 95% CIs
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(exp(fixedres[ , grep(pattern="beta0.", x=names(fixedres))]),
                outline=FALSE, col='gray87', ylim=c(0,300),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    exp(fixedres[ , grep(pattern="beta0.", x=names(fixedres))]),
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,300),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = FALSE)
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  #mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic(omega)[g]), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="sw", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )  
  
# Estimation accuracy for omega from random effects model
# Compare estimated omega to true values using 
# a boxplot that shows 95% CIs
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(exp(ranres[ , grep(pattern="beta0.", x=names(ranres))]),
                outline=FALSE, col='gray87', ylim=c(0,300),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    exp(ranres[ , grep(pattern="beta0.", x=names(ranres))]),
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,300),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = FALSE)
  axis(side=2, las=2, labels=FALSE)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  #mtext(text='Population number', side=1, line=3.5)
  #mtext(text=expression(italic(omega)[g]), side=2, line=3.5)
  points(1:10,
         ranres[1 , grep(pattern="sw", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )  
    
# Estimation accuracy for t0 from fixed effects model
# Compare estimated t0 to true values using 
# a boxplot that shows 95% CIs
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(fixedres[ , grep(pattern="to.", x=names(fixedres))],
                outline=FALSE, col='gray87', ylim=c(-2,0),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    fixedres[ , grep(pattern="to.", x=names(fixedres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(-2,0),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  #mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression('t'['0'['g']]), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="st0", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )
  
# Estimation accuracy for t0 from random effects model
# Compare estimated t0 to true values using 
# a boxplot that shows 95% CIs
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ranres[ , grep(pattern="to.", x=names(ranres))],
                outline=FALSE, col='gray87', ylim=c(-2,0),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ranres[ , grep(pattern="to.", x=names(ranres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(-2,0),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=1,
      medcol='gray40', medlwd=1)
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2, labels=FALSE)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5, adj=-1.25)

  #mtext(text=expression('t'['0'['g']]), side=2, line=3.5)
  points(1:10,
         ranres[1 , grep(pattern="st0", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )
  dev.off()