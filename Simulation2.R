library(lars)
library(grpreg)
library(KernSmooth)
library(DNAcopy)    #download from https://bioconductor.org/packages/release/bioc/html/DNAcopy.html

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

localDiagnostic<-function (y, h) {    #see original SaRa algorithm 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

localMax<-function (y, span = 5) {   #see original SaRa algorithm
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

refine<-function (a, error=2){   #delete adjacent estimators
  x=a[1]
  if(length(a)>1){
    for(i in 2:length(a)){
      if(min(abs(x-a[i]))>error)
        x=c(x,a[i])
    }
  }
  return(x)
}

######simulation 2
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
#######Case 2, d=1
set.seed(23)
TP1=matrix(NA,nr=100,nc=14)                #true positive
TP2=matrix(NA,nr=100,nc=14)
TP3=matrix(NA,nr=100,nc=14)
TP4=matrix(NA,nr=100,nc=14)
TP5=matrix(NA,nr=100,nc=14)
FP1=matrix(NA,nr=100,nc=14)                #false positive
FP2=matrix(NA,nr=100,nc=14)
FP3=matrix(NA,nr=100,nc=14)
FP4=matrix(NA,nr=100,nc=14)
FP5=matrix(NA,nr=100,nc=14)
Haus1=matrix(NA,nr=100,nc=14)              #hausdorff distance
Haus2=matrix(NA,nr=100,nc=14)
Haus3=matrix(NA,nr=100,nc=14)
Haus4=matrix(NA,nr=100,nc=14)
Haus5=matrix(NA,nr=100,nc=14)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=arima.sim(n=1100, list(ar = 0.2), sd = 0.2)[101:1100]    
  y=signal+e  
  ######CBS
  CBS=DNAcopy::segment(CNA(y, rep(1,n), 1:n))
  num1=length(CBS$output[,4])-1
  estimate=refine(CBS$output[1:num1,4])
  TP1[i,]=TruePositive(estimate, tau)
  FP1[i,]=FalsePositive(estimate, tau)
  Haus1[i,]=hausdorff(estimate, tau)
  ######SaRa
  SaRa=abs(localDiagnostic(y, h=10))
  candidate=localMax(SaRa)
  index=candidate[order(SaRa[candidate],decreasing=TRUE)]     #ranking
  for (j in 1:14){
    estimate=index[1:j]
    TP2[i,length(estimate)]=TruePositive(estimate, tau)
    FP2[i,length(estimate)]=FalsePositive(estimate, tau)
    Haus2[i,length(estimate)]=hausdorff(estimate, tau)
  }
  ######standard fused-Lasso
  model1=lars(X, y=y, type="lasso", normalize=FALSE, intercept=TRUE, trace=FALSE, max.steps=30)
  estimate=refine(unlist(model1$actions))
  for(j in 1:14){
    TP3[i,j]=TruePositive(estimate[1:j], tau)
    FP3[i,j]=FalsePositive(estimate[1:j], tau)
    Haus3[i,j]=hausdorff(estimate[1:j], tau)
  }
  ######proposed I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=30)
  estimate=refine(unlist(model2$actions))
  for(j in 1:14){
    TP4[i,j]=TruePositive(estimate[1:j], tau)
    FP4[i,j]=FalsePositive(estimate[1:j], tau)
    Haus4[i,j]=hausdorff(estimate[1:j], tau)
  }
  ######proposed II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=30)
  estimate=refine(unlist(model3$actions))
  for(j in 1:14){
    TP5[i,j]=TruePositive(estimate[1:j], tau)
    FP5[i,j]=FalsePositive(estimate[1:j], tau)
    Haus5[i,j]=hausdorff(estimate[1:j], tau)
  }
}
######show Figure 
par(mfrow=c(2,3))
plot(x=1:14,y=apply(TP2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,11), xlab = "Number of selected change points",ylab ="True positives",main="TP, d=1")
lines(x=1:14,y=apply(TP3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(TP4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(TP5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(TP1,2,mean,na.rm=T),type="l",lty=2)
#legend("bottomright", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")
plot(x=1:14,y=apply(FP2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,11), xlab = "Number of selected change points",ylab ="False positives",main="FP, d=1")
lines(x=1:14,y=apply(FP3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(FP4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(FP5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(FP1,2,mean,na.rm=T),type="l",lty=2)
#legend("topleft", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")
plot(x=1:14,y=apply(Haus2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff, d=1")
lines(x=1:14,y=apply(Haus3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(Haus4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(Haus5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(Haus1,2,mean,na.rm=T),type="l",lty=2)
legend("topright", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")

#######Case 2, d=4
d=4
group=NULL
for (j in 1:n) group<-c(group,rep(j-1,d))
XXX=kronecker(x1,diag(1, d))
XX1=kronecker(cbind(x1[,1],X1), diag(1, d))
XX2=kronecker(cbind(x1[,1],X2), diag(1, d))
TP1=matrix(NA,nr=100,nc=14)                #true positive
TP2=matrix(NA,nr=100,nc=14)
TP3=matrix(NA,nr=100,nc=14)
TP4=matrix(NA,nr=100,nc=14)
TP5=matrix(NA,nr=100,nc=14)
FP1=matrix(NA,nr=100,nc=14)                #false positive
FP2=matrix(NA,nr=100,nc=14)
FP3=matrix(NA,nr=100,nc=14)
FP4=matrix(NA,nr=100,nc=14)
FP5=matrix(NA,nr=100,nc=14)
Haus1=matrix(NA,nr=100,nc=14)              #hausdorff distance
Haus2=matrix(NA,nr=100,nc=14)
Haus3=matrix(NA,nr=100,nc=14)
Haus4=matrix(NA,nr=100,nc=14)
Haus5=matrix(NA,nr=100,nc=14)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE) 
  beta=matrix(sample(c(-1,-0.5,0.5,1), size=J*d, replace = TRUE), nrow=d)
  signal=matrix(rep(0, n*d), nrow=d)
  y=matrix(rep(0, n*d), nrow=d)
  for(k in 1:d){
    signal[k, ]=0
    for(j in 1:J) {signal[k, ]<-signal[k, ]+beta[k,j]*(x>tau[j])}
    y[k, ]=signal[k, ]+arima.sim(n=1100, list(ar = 0.2), sd = 0.2)[101:1100]
  }
  ######CBS
  CBS=DNAcopy::segment(CNA(t(y), chrom=rep(1,n), maploc=1:n))
  candidate=sort(unique(CBS$output[,4]))
  estimate=refine(candidate[-length(candidate)])
  TP1[i,]=TruePositive(estimate, tau)
  FP1[i,]=FalsePositive(estimate, tau)
  Haus1[i,]=hausdorff(estimate, tau)
  ######SaRa
  SaRa=localDiagnostic(y[1,], h=10)^2+localDiagnostic(y[2,], h=10)^2+localDiagnostic(y[3,], h=10)^2+localDiagnostic(y[4,], h=10)^2
  candidate=localMax(SaRa,span=10)
  index=candidate[order(SaRa[candidate],decreasing=TRUE)]     #ranking
  for (j in 1:14){
    estimate=index[1:j]
    TP2[i,length(estimate)]=TruePositive(estimate, tau)
    FP2[i,length(estimate)]=FalsePositive(estimate, tau)
    Haus2[i,length(estimate)]=hausdorff(estimate, tau)
  }
  ######standard fused-Lasso
  yy=as.vector(y)                 
  gLasso=grpreg(X=XXX[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grLasso",family="gaussian", gmax=20)   
  for(j in 2:dim(gLasso$beta)[2]){
    candidate=as.vector((which(gLasso$beta[,j]!=0)[d*(1:(length(which(gLasso$beta[,j]!=0))/d))]/d-1)[-1])
    estimate=refine(candidate)
    if(length(estimate)<=14){
      TP3[i,length(estimate)]=TruePositive(estimate, tau)
      FP3[i,length(estimate)]=FalsePositive(estimate, tau)
      Haus3[i,length(estimate)]=hausdorff(estimate, tau) 
    }
  }
  ######Proposed I
  yy=as.vector(matrix(c(y[1,]-locpoly(x,y[1,],kernel="epanech",bandwidth=h,gridsize=n)$y, y[2,]-locpoly(x,y[2,],kernel="epanech",bandwidth=h,gridsize=n)$y, 
                        y[3,]-locpoly(x,y[3,],kernel="epanech",bandwidth=h,gridsize=n)$y, y[4,]-locpoly(x,y[4,],kernel="epanech",bandwidth=h,gridsize=n)$y), nrow=d, byrow=T))
  model1=grpreg(X=XX1[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grLasso",family="gaussian", gmax=20) 
  for(j in 2:dim(model1$beta)[2]){
    candidate=as.vector((which(model1$beta[,j]!=0)[d*(1:(length(which(model1$beta[,j]!=0))/d))]/d-1)[-1])
    estimate=refine(candidate)
    if(length(estimate)<=14){
      TP4[i,length(estimate)]=TruePositive(estimate, tau)
      FP4[i,length(estimate)]=FalsePositive(estimate, tau)
      Haus4[i,length(estimate)]=hausdorff(estimate, tau) 
    }
  }
  ######Proposed II
  yy=as.vector(matrix(c(eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%(y[1,]-mean(y[1,])), eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%(y[2,]-mean(y[2,])),
                        eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%(y[3,]-mean(y[3,])), eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%(y[4,]-mean(y[4,]))), nrow=d, byrow=T))
  model2=grpreg(X=XX2[,2:length(yy)], y=yy, group=group[2:length(yy)], penalty="grLasso",family="gaussian", gmax=20) 
  for(j in 2:dim(model2$beta)[2]){
    candidate=as.vector((which(model2$beta[,j]!=0)[d*(1:(length(which(model2$beta[,j]!=0))/d))]/d-1)[-1])
    estimate=refine(candidate)
    if(length(estimate)<=14){
      TP5[i,length(estimate)]=TruePositive(estimate, tau)
      FP5[i,length(estimate)]=FalsePositive(estimate, tau)
      Haus5[i,length(estimate)]=hausdorff(estimate, tau) 
    }
  }
}
######show Figure
plot(x=1:14,y=apply(TP2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,11), xlab = "Number of selected change points",ylab ="True positives",main="TP, d=4")
lines(x=1:14,y=apply(TP3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(TP4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(TP5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(TP1,2,mean,na.rm=T),type="l",lty=2)
#legend("bottomright", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")
plot(x=1:14,y=apply(FP2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,11), xlab = "Number of selected change points",ylab ="False positives",main="FP, d=4")
lines(x=1:14,y=apply(FP3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(FP4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(FP5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(FP1,2,mean,na.rm=T),type="l",lty=2)
#legend("topleft", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")
plot(x=1:14,y=apply(Haus2,2,mean,na.rm=T),pch=1,col=2, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff, d=4")
lines(x=1:14,y=apply(Haus3,2,mean,na.rm=T),pch=2,col=3, type="b")
lines(x=1:14,y=apply(Haus4,2,mean,na.rm=T),pch=3,col=4, type="b")
lines(x=1:14,y=apply(Haus5,2,mean,na.rm=T),pch=4,col=5, type="b")
lines(x=1:14,y=apply(Haus1,2,mean,na.rm=T),type="l",lty=2)
legend("topright", legend=c("SaRa","fused-Lasso","Proposed I","Proposed II"), pch=c(1,2,3,4), col=c(2,3,4,5), bty="n")
