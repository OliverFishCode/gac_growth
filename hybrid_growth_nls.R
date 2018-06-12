setwd("c:/Users/olive/Google Drive/growth lecture")#sets working directory
Set1 = data.frame(read.csv(file="C:/Users/olive/Google Drive/growth lecture/set1.csv"))#calls dataset
Set2 = data.frame(read.csv(file="C:/Users/olive/Google Drive/growth lecture/set2.csv"))#calls dataset
Set3 = rbind(Set1,Set2)
library(nlstools)

#############Von Bertalanffy Growth Model##############
vonb = nls((tl ~ linf * (1- exp(-k*(age - t0)))),start=list(linf=400, k=0.73,t0=0), data=Set3) #model,intial values, data set... fits growth model
vbp = data.frame(predict(vonb, interval = confidence)) #provide model predictions 
summary(vonb) #produces parameter estmates, std. error, t value and associated p-values
nlstools::confint2(vonb) #produces 95% confidence interval for parameter estimates

#############Gallucci and Quinn formulation of Von Bertalanffy Growth Model ##############
############# Linf is replaced by w/k, where w= linf*k##################################
GallQuinn = nls((tl ~ w/k* (1- exp(-k*(age - t0)))),start=list(w=292, k=0.73,t0=0), data=Set3) 
Gp = data.frame(predict(GallQuinn))
summary(GallQuinn) 
nlstools::confint2(GallQuinn)

#############Ogle and Iserman Von Bertalanffy Growth Model- mean length at specified age##############
Ogleisermann = nls((tl ~ lr + (linf - lr) * (1- exp(-k*(age - 2)))),start=list(linf=400, k=0.73,lr=200), data=Set3) 
oip = data.frame(predict(Ogleisermann)) 
summary(Ogleisermann ) 
nlstools::confint2(Ogleisermann) 

#############Ogle and Isermann Von Bertalanffy Growth Model- age at specifed length##############
Ogleisermann2  = nls((tl ~ 213.95574 + (linf - 213.95574) * (1- exp(-k*(age - tr)))),start=list(linf=400, k=0.73,tr=1), data=Set3) 
oi2p = data.frame(predict(Ogleisermann2)) 
summary(Ogleisermann2) 
nlstools::confint2(Ogleisermann2)

#############Oi_gac hybrid lr estimation##############
hybrid_lr = nls((tl ~ lr + ((w/k) - lr) * (1- exp(-k*(age - 0)))),start=list(w=100, k=0.73,lr=77), data=Set3) 
hyblr_p = data.frame(predict(hybrid_lr)) 
summary(hybrid_lr) 
nlstools::confint2(hybrid_lr) 

#############Oi_gac hybrid tr estimation##############
hybrid_tr = nls((tl ~77.59 + ((w/k) - 77.59) * (1- exp(-k*(age - tr)))),start=list(w=292, k=0.73, tr=0), data=Set3) 
hybtr_p = data.frame(predict(hybrid_tr)) 
summary(hybrid_tr) 
nlstools::confint2(hybrid_tr)


#############Combine Predictions##############
age = subset(Set3,select= -c(tl,id))
allp = data.frame(vbp,Gp,oip,oi2p,hyblr_p,hybtr_p,age)