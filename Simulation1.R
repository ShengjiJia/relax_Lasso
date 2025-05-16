library(lars)
library(KernSmooth)

######define some functions
TruePositive<-function(estimate, true, error=2){         
  #calculate the number of true positives
  Num=0
  if(length(estimate)>0){
    for(l in 1:length(true)){
      if(min(abs(estimate-true[l]))<=error)
        Num=Num+1
    }
  }
  return(Num)
}

FalsePositive<-function(estimate, true, error=2){         
  #calculate the number of flase positives
  Num=0
  if(length(estimate)>0){
    for(l in 1:length(estimate)){
      if(min(abs(estimate[l]-true))>error)
        Num=Num+1
    }
  }
  return(Num)
}

hausdorff<-function(x, y){
  if (length(y)*length(x)==0){  
    return(NA)
  }   
  else{
    A=rep(0, length(y))
    for(i in 1:length(y)){
      A[i]=min(abs(y[i]-x))  
    }
    B=rep(0, length(x))
    for(i in 1:length(x)){
      B[i]=min(abs(x[i]-y))
    }
    return(max(max(A),max(B)))
  }
}
######simulation 1
n=1000
x=1:n
J=10
h=20
lambda=0.5
x1=matrix(1,nrow=n,ncol=n)      
for(i in 1:(n-1)) {x1[i,(i+1):n]=0}
X=x1[ ,2:n]                               #matrix X
X1=X                                       #matrix (I-S)X
for(i in 1:(n-1)) {X1[,i]=X[,i]-locpoly(x,X[,i],kernel="epanech",bandwidth=h,gridsize=n)$y}
XX=X                                       #centered X
for(i in 1:(n-1)) {XX[,i]=X[,i]-mean(X[,i])}
R=XX%*%solve(t(XX)%*%XX+n*lambda*diag(1,n-1))%*%t(XX)
eig=eigen(diag(1,n)-R)
X2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%XX       #matrix (I-R)^{1/2}(I-11/n)X
TP0=matrix(NA,nr=100,nc=15)                #true positive
TP1=matrix(NA,nr=100,nc=15)
tp1=matrix(NA,nr=100,nc=15)
FP0=matrix(NA,nr=100,nc=15)                #false positive
FP1=matrix(NA,nr=100,nc=15)
fp1=matrix(NA,nr=100,nc=15)
Haus0=matrix(NA,nr=100,nc=15)              #hausdorff distance
Haus1=matrix(NA,nr=100,nc=15)
haus1=matrix(NA,nr=100,nc=15)
set.seed(18)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=rnorm(n,mean=0,sd=0.2)      
  y=signal+e  
  ######standard fused-Lasso
  model1=lars(X, y=y, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate1=unlist(model1$actions)[1:j]
    TP0[i,length(estimate1)]=TruePositive(estimate1, tau)
    FP0[i,length(estimate1)]=FalsePositive(estimate1, tau)
    Haus0[i,length(estimate1)]=hausdorff(estimate1, tau)
  }
  ######proposed I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP1[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP1[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus1[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######proposed II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    tp1[i,length(estimate3)]=TruePositive(estimate3, tau)
    fp1[i,length(estimate3)]=FalsePositive(estimate3, tau)
    haus1[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
}

######show Figure 
par(mfrow=c(1,3))
plot(x=1:15,y=apply(TP0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(tp1,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("fused-Lasso","Proposed I","Proposed II"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(fp1,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("fused-Lasso","Proposed I","Proposed II"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(haus1,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("fused-Lasso","Proposed I","Proposed II"), pch=c(2,1,4), col=c(3,2,4), bty="n")

tau
estimate1
estimate2
estimate3

######compare different tuning parameters
h=40
lambda=1
X1=X                                       #matrix (I-S)X
for(i in 1:(n-1)) {X1[,i]=X[,i]-locpoly(x,X[,i],kernel="epanech",bandwidth=h,gridsize=n)$y}
XX=X                                       #centered X
for(i in 1:(n-1)) {XX[,i]=X[,i]-mean(X[,i])}
R=XX%*%solve(t(XX)%*%XX+n*lambda*diag(1,n-1))%*%t(XX)
eig=eigen(diag(1,n)-R)
X2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%XX       #matrix (I-R)^{1/2}(I-11/n)X
TP2=matrix(NA,nr=100,nc=15)
tp2=matrix(NA,nr=100,nc=15)
FP2=matrix(NA,nr=100,nc=15)
fp2=matrix(NA,nr=100,nc=15)
Haus2=matrix(NA,nr=100,nc=15)
haus2=matrix(NA,nr=100,nc=15)
set.seed(18)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=rnorm(n,mean=0,sd=0.2)      
  y=signal+e  
  ######proposed I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP2[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP2[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus2[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######proposed II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    tp2[i,length(estimate3)]=TruePositive(estimate3, tau)
    fp2[i,length(estimate3)]=FalsePositive(estimate3, tau)
    haus2[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
}

h=10
lambda=0.25
X1=X                                       #matrix (I-S)X
for(i in 1:(n-1)) {X1[,i]=X[,i]-locpoly(x,X[,i],kernel="epanech",bandwidth=h,gridsize=n)$y}
XX=X                                       #centered X
for(i in 1:(n-1)) {XX[,i]=X[,i]-mean(X[,i])}
R=XX%*%solve(t(XX)%*%XX+n*lambda*diag(1,n-1))%*%t(XX)
eig=eigen(diag(1,n)-R)
X2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%XX       #matrix (I-R)^{1/2}(I-11/n)X
TP3=matrix(NA,nr=100,nc=15)
tp3=matrix(NA,nr=100,nc=15)
FP3=matrix(NA,nr=100,nc=15)
fp3=matrix(NA,nr=100,nc=15)
Haus3=matrix(NA,nr=100,nc=15)
haus3=matrix(NA,nr=100,nc=15)
set.seed(18)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=rnorm(n,mean=0,sd=0.2)      
  y=signal+e  
  ######proposed I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP3[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP3[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus3[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######proposed II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    tp3[i,length(estimate3)]=TruePositive(estimate3, tau)
    fp3[i,length(estimate3)]=FalsePositive(estimate3, tau)
    haus3[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
}

######show Figure 
par(mfrow=c(2,3))
plot(x=1:15,y=apply(TP1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(tp1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(tp2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(tp3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("lambda2=0.5","lambda2=1","lambda2=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(fp1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(fp2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(fp3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("lambda2=0.5","lambda2=1","lambda2=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(haus1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(haus2,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(haus3,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("lambda2=0.5","lambda2=1","lambda2=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")

