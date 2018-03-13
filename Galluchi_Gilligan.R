library(FSA)
library(nlstools)
library(dplyr)
fish <-read.csv('set2.csv')
fish2<-read.csv('set1.csv')

######
# Estimated Parameters using Gallucci and Quinn for set #2
headtail(fish)
(st.values<-vbStarts(tl~age,data=fish))
vb1<-vbFuns("GallucciQuinn")
fitgq<-nls(tl~vb1(age,Linf,K,t0),start=st.values, data=fish)
summary(fitgq,correlation=TRUE)
coef(fitgq)
confint(fitgq)
#Plot
xf<-seq(min(fish$age),max(fish$age),length.out=199)
pf<-vb1(xf,omega=coef(fitgq))
xlmts<-range(c(xf))
ylmts<-range(c(fish$tl))
plot(tl~age,data=fish,xlab="Age",ylab="Total Length",xlim=xlmts,ylim=ylmts,col="White")
points(tl~age,data=fish,pch=1,col=rgb(0,0,0,1/2),cex=0.8)
lines(pf~xf,lwd=2,lty="solid")
# Predictions of lengths
nd<-data.frame(age=c(1,2,3,4,5,6))
predict(fitgq,nd)
#####
# Estimated Parameters using Gallucci and Quinn for set #1
headtail(fish2)
(st.values2<-vbStarts(tl~age,data=fish2))
vb2<-vbFuns("GallucciQuinn")
fitgq2<-nls(tl~vb2(age,Linf,K,t0),start=st.values2, data=fish2)
summary(fitgq2,correlation=TRUE)
coef(fitgq2)
confint(fitgq2)
#Plot
xf2<-seq(min(fish2$age),max(fish2$age),length.out=199)
pf2<-vb2(xf2,omega=coef(fitgq2))
xlmts2<-range(c(xf2))
ylmts2<-range(c(fish2$tl))
plot(tl~age,data=fish2,xlab="Age",ylab="Total Length",xlim=xlmts2,ylim=ylmts2,col="White")
points(tl~age,data=fish2,pch=2,col=rgb(0,0,0,1/2),cex=0.8)
lines(pf2~xf2,lwd=2,lty="dashed")
# Predictions of lengths
nd2<-data.frame(age=c(1,2,3,4,5,6))
predict(fitgq2,nd2)

#plot containing both
plot(tl~age,data=fish2,xlab="Age",ylab="Total Length",xlim=xlmts2,ylim=ylmts2,col="White")
points(tl~age,data=fish,pch=20,col=("blue"),cex=0.8)
lines(pf~xf,lwd=2,col="blue",lty="solid")
points(tl~age,data=fish2,pch=21,col="red",cex=0.8)
lines(pf2~xf2,lwd=2,col="red",lty="dashed")
