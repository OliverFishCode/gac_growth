# Base model -----
load('result.rda')
head(res)
nrow(res)
# Estimation accuracy for k
par(mar=c(5,5,1,1))
hist(res$K, col='gray87', xlab=expression(italic('k')),
     ylab = 'Frequency', yaxt='n', xaxt='n',
     xlim = c(0, .5), main='')
abline(v=mean(res$k), col='blue', lty=1, lwd=3)
axis(side=1, pos=0)    
axis(side=2, pos=0, las=TRUE)

# Estimation accuracy for w
par(mar=c(5,5,1,1))
hist(exp(res$beta0), col='gray87', 
     xlab=expression(omega),
     ylab = 'Frequency', yaxt='n', xaxt='n',
     xlim = c(80, 160), main='')
abline(v=mean(res$k)*mean(res$linf), col='blue', lty=1, lwd=3)
axis(side=1, pos=0)    
axis(side=2, pos=80, las=TRUE)

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


