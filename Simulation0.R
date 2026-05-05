library(lars)
library(glmnet)
library(grpreg)
library(KernSmooth)

######simulation in introduction
set.seed(6)
n=1000
x=1:n
x1=matrix(1,nrow=n,ncol=n)     #matrix X
for(i in 1:(n-1)) x1[i,(i+1):n]=0
J=6
tau=c(50,300,450,500,700,800)
beta=c(-2,1.5,-2.5,1,2,-1.5)
signal=2
for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
e=rnorm(n,mean=0,sd=0.5)      
y=signal+e  
######lars
model1=lars(x1[1:n,2:n], y=y, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=16)
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
######relax fused-Lasso I
h=20
x10=x1
for(i in 2:n) {x10[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }  
y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
model3=lars(x10[1:n,2:n], y=y0, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=10)
estimate3=unlist(model3$actions)
model3=lars(x10[1:n,2:n], y=y0, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=9)
estimate=unlist(model3$actions)
Estimate3=c(0,sort(estimate),n)
Fit3=NULL
for(i in 1:(length(estimate)+1)){
  m=mean(y[(Estimate3[i]+1):Estimate3[i+1]])
  Fit3=c(Fit3,rep(m, Estimate3[i+1]-Estimate3[i]))
}
######relax fused-SCAD I
model4=grpreg(X=x10[1:n,2:n], y=y0, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=8) 
estimate4=as.vector((which(model4$beta[2:n,1]!=0)))
for(j in 2:dim(model4$beta)[2]){
  s=setdiff(as.vector((which(model4$beta[2:n,j]!=0))), estimate4)
  if(length(s)>0) estimate4=c(estimate4, s)
}
######show Figure1
par(mfrow=c(4,1), mar=c(3,4,2,1))
plot(x, y, ylim=c(-1.5,2.5), xlab="", ylab="signal", main="true signal and data", pch=20, col=8)
lines(x, signal, type="l", col=2, lwd=2)
abline(v=tau, lty=2)
plot(x, Fit1, type="l", ylim=c(-1.5,2.5), xlab="", ylab="signal", main="fused-Lasso", col=2, lwd=2)
abline(v=tau, lty=2)
plot(x, Fit2, type="l", ylim=c(-1.5,2.5), xlab="", ylab="signal", main="elastic net", col=2, lwd=2)
abline(v=tau, lty=2)
plot(x, Fit3, type="l", ylim=c(-1.5,2.5), xlab="locations", ylab="signal", main="relaxed fused-Lasso", col=2, lwd=2)
abline(v=tau, lty=2)

estimate1
estimate3
estimate4
