library(lars)
library(glmnet)
library(KernSmooth)

######simulation in introduction
set.seed(18)
n=1000
x=1:n
x1=matrix(1,nrow=n,ncol=n)     #matrix X
for(i in 1:(n-1)) x1[i,(i+1):n]=0
tau=c(150,350,600,850)
beta=c(1,-1,-1,1)
signal=0
J=4
for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
e=rnorm(n,mean=0,sd=0.5)      
y=signal+e  
######lars
model1=lars(x1[1:n,2:n], y=y, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=10)
estimate1=unlist(model1$actions)
Estimate1=c(0,sort(estimate1),n)
Fit1=NULL
for(i in 1:(length(estimate1)+1)){
  m=mean(y[(Estimate1[i]+1):Estimate1[i+1]])
  Fit1=c(Fit1,rep(m, Estimate1[i+1]-Estimate1[i]))
}
######elestic net
model2=glmnet(x1[1:n,2:n], y=y, family="gaussian", alpha=0.5, standardize=F, lambda=0.02)
estimate2=which(model2$beta!=0)
Estimate2=c(0,estimate2,n)
Fit2=NULL
for(i in 1:(length(estimate2)+1)){
  m=mean(y[(Estimate2[i]+1):Estimate2[i+1]])
  Fit2=c(Fit2,rep(m, Estimate2[i+1]-Estimate2[i]))
}
######proposed
h=20
x10=x1
for(i in 2:n) {x10[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }  
y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
#model3=glmnet(x10[1:n,2:n], y=y0, family="gaussian", alpha=1, standardize=F, lambda=0.003)
#estimate3=which(model3$beta!=0)
model3=lars(x10[1:n,2:n], y=y0, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=5)
estimate3=unlist(model3$actions)
Estimate3=c(0,sort(estimate3),n)
Fit3=NULL
for(i in 1:(length(estimate3)+1)){
  m=mean(y[(Estimate3[i]+1):Estimate3[i+1]])
  Fit3=c(Fit3,rep(m, Estimate3[i+1]-Estimate3[i]))
}
######show Figure1
par(mfrow=c(4,1), mar=c(3,4,2,1))
plot(x, y, ylim=c(-2,2), xlab="", ylab="signal", main="true signal and data", pch=20, col=8)
lines(x, signal, type="l", col=2, lwd=2)
abline(v=c(150,350,600,850), lty=2)
plot(x, Fit1, type="l", ylim=c(-1.5,1.5), xlab="", ylab="signal", main="fused-Lasso", col=2, lwd=2)
abline(v=c(150,350,600,850), lty=2)
plot(x, Fit2, type="l", ylim=c(-1.5,1.5), xlab="", ylab="signal", main="elastic net", col=2, lwd=2)
abline(v=c(150,350,600,850), lty=2)
plot(x, Fit3, type="l", ylim=c(-1.5,1.5), xlab="locations", ylab="signal", main="relaxed fused-Lasso", col=2, lwd=2)
abline(v=c(150,350,600,850), lty=2)

estimate1
estimate3
