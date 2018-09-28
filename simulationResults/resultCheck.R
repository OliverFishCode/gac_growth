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
       xlim = c(0, 1), main='')
  abline(v=mean(covres$sk), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=0, las=TRUE)

# Estimation accuracy for betaT
  par(mar=c(5,5,1,1))
  hist(covres$betaT, col='gray87', 
       xlab=expression(beta['T']),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(-.6, .6), main='')
  abline(v=mean(covres$sbetaT), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=-.6, las=TRUE)

# Estimation accuracy for w
# Calculate w as function of linear predictor
  ew = exp(covres$beta0 + covres$betaT*covres$temp)
# Make plot
  par(mar=c(5,5,1,1))
  hist(ew, col='gray87', 
       xlab=expression(omega),
       ylab = 'Frequency', yaxt='n', xaxt='n',
       xlim = c(0, 300), main='')
  abline(v=mean(covres$sw), col='blue', lty=1, lwd=3)
  axis(side=1, pos=0)    
  axis(side=2, pos=0, las=TRUE)

# Fixed-effects model -----
  load('fixedresult.rda')
  head(fixedres)
  nrow(fixedres)  
  
# Estimation accuracy for k
# Compare estimated k to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(fixedres[ , grep(pattern="K.", x=names(fixedres))],
                outline=FALSE, col='gray87', ylim=c(0,1),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    fixedres[ , grep(pattern="K.", x=names(fixedres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,1),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic('k'[g])), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="sk", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )

# Estimation accuracy for omega
# Compare estimated omega to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
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
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic(omega[g])), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="sw", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )  
  
# Estimation accuracy for t0
# Compare estimated t0 to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
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
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic('t')['0,'][italic(' g')]), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="st0", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )
  
# . Fig. 2 ----
  # Set up graphics device
  tiff(
    filename = "Figure2.tif",
    width = 970,
    height = 1280,
    pointsize = 8,
    res = 350
  )

 # Set up plotting window
    par(mfrow=c(3, 1), oma=c(4,1.1,1,1), mar=c(1,4,1,1))
 # Estimation accuracy for k
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
  #mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression('K'[g]), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="sk", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )

# Estimation accuracy for omega
# Compare estimated omega to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(1,4,1,1))
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
  
# Estimation accuracy for t0
# Compare estimated t0 to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(1,4,1,1))
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
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression('t'['0'['g']]), side=2, line=3.5)
  points(1:10,
         fixedres[1 , grep(pattern="st0", x=names(fixedres))],
         pch=21, col='white', bg='gray40'
    )
  dev.off()
  
# Random-effects model -----
  load('ranresult.rda')
  head(ranres)
  nrow(ranres)  
  
# Estimation accuracy for k
# Compare estimated k to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ranres[ , grep(pattern="K.", x=names(ranres))],
                outline=FALSE, col='gray87', ylim=c(0,1),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ranres[ , grep(pattern="K.", x=names(ranres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(0,1),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic('k'[g])), side=2, line=3.5)
  points(1:10,
         ranres[1 , grep(pattern="sk", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )

# Estimation accuracy for omega
# Compare estimated omega to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
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
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic(omega[g])), side=2, line=3.5)
  points(1:10,
         ranres[1 , grep(pattern="sk", x=names(ranres))]*
           ranres[1 , grep(pattern="slinf", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )  
  
# Estimation accuracy for t0
# Compare estimated t0 to true values using 
# a boxplot that shows 95% CIs
  # Set margins
  par(mar=c(5,5,1,1))
  # Set up the boxplot, will write over this with new whiskers
  bb <- boxplot(ranres[ , grep(pattern="to.", x=names(ranres))],
                outline=FALSE, col='gray87', ylim=c(-10,0),
                col.axis='white', notch=FALSE, plot=FALSE)
  # Replace stats for whiskers with 95% CI
  bb$stats[c(1,5), ] <- apply(
    ranres[ , grep(pattern="to.", x=names(ranres))],
    2,
    quantile,
    probs=c(.025, 0.975), na.rm = TRUE
    )
  # Re-plot
  bxp(bb, outline=FALSE, ylim=c(-10,0),
      boxfill='gray87', col.axis='white',
      staplewex=0, whisklty=1, whiskcol='gray40',
      whisklwd=2, boxcol='gray40', boxlwd=2,
      medcol='gray40')
  # Add x (side=1) and y (side=2) axes
  axis(side=1, at=seq(1,10,1), labels = seq(1,10,1))
  axis(side=2, las=2)
  # Close it up with a box so axis styles match
  box()
  # Add axis labels
  mtext(text='Population number', side=1, line=3.5)
  mtext(text=expression(italic('t')['0,'][italic(' g')]), side=2, line=3.5)
  points(1:10,
         ranres[1 , grep(pattern="st0", x=names(ranres))],
         pch=21, col='white', bg='gray40'
    )
  