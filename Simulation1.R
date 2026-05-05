library(lars)
library(KernSmooth)
library(grpreg)

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
TP2=matrix(NA,nr=100,nc=15)
TP3=matrix(NA,nr=100,nc=15)
TP4=matrix(NA,nr=100,nc=15)
FP0=matrix(NA,nr=100,nc=15)                #false positive
FP1=matrix(NA,nr=100,nc=15)
FP2=matrix(NA,nr=100,nc=15)
FP3=matrix(NA,nr=100,nc=15)
FP4=matrix(NA,nr=100,nc=15)
Haus0=matrix(NA,nr=100,nc=15)              #hausdorff distance
Haus1=matrix(NA,nr=100,nc=15)
Haus2=matrix(NA,nr=100,nc=15)
Haus3=matrix(NA,nr=100,nc=15)
Haus4=matrix(NA,nr=100,nc=15)
set.seed(34)
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
  ######relax fused-Lasso I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP1[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP1[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus1[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######relax fused-Lasso II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    TP2[i,length(estimate3)]=TruePositive(estimate3, tau)
    FP2[i,length(estimate3)]=FalsePositive(estimate3, tau)
    Haus2[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
  #######relax fused-SCAD I
  model4=grpreg(X=X1, y=y1, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model4$beta)[2]){
    estimate4=as.vector((which(model4$beta[2:n,j]!=0)))
    if(length(estimate4)<16){
      TP3[i,length(estimate4)]=TruePositive(estimate4, tau)
      FP3[i,length(estimate4)]=FalsePositive(estimate4, tau)
      Haus3[i,length(estimate4)]=hausdorff(estimate4, tau)
    }
  }
  ########relax fused-SCAD II
  model5=grpreg(X=X2, y=y2, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model5$beta)[2]){
    estimate5=as.vector((which(model5$beta[2:n,j]!=0)))
    if(length(estimate5)<16){
      TP4[i,length(estimate5)]=TruePositive(estimate5, tau)
      FP4[i,length(estimate5)]=FalsePositive(estimate5, tau)
      Haus4[i,length(estimate5)]=hausdorff(estimate5, tau)
    }
  }
}
estimate4=as.vector((which(model4$beta[2:n,1]!=0)))
for(j in 2:dim(model4$beta)[2]){
  s=setdiff(as.vector((which(model4$beta[2:n,j]!=0))), estimate4)
  if(length(s)>0) estimate4=c(estimate4, s)
}
estimate5=as.vector((which(model5$beta[2:n,1]!=0)))
for(j in 2:dim(model5$beta)[2]){
  s=setdiff(as.vector((which(model5$beta[2:n,j]!=0))), estimate5)
  if(length(s)>0) estimate5=c(estimate5, s)
}

######show Figure 
par(mfrow=c(1,3))
plot(x=1:15,y=apply(TP0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP2,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:15,y=apply(TP3,2,mean,na.rm=T),pch=5, col=7, type="b")
lines(x=1:15,y=apply(TP4,2,mean,na.rm=T),pch=6, col=6, type="b")
#legend("bottomright", legend=c("fused-Lasso","relex fused-Lasso I","relex fused-Lasso II", "relex fused-SCAD I", "relex fused-SCAD II"), pch=c(2,1,4,5,6), col=c(3,2,4,7,6), bty="n")
plot(x=1:15,y=apply(FP0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP2,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:15,y=apply(FP3,2,mean,na.rm=T),pch=5, col=7, type="b")
lines(x=1:15,y=apply(FP4,2,mean,na.rm=T),pch=6, col=6, type="b")
legend("topleft", legend=c("fused-Lasso","relex fused-Lasso I","relex fused-Lasso II", "relex fused-SCAD I", "relex fused-SCAD II"), pch=c(2,1,4,5,6), col=c(3,2,4,7,6), bty="n")
plot(x=1:15,y=apply(Haus0,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus1,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus2,2,mean,na.rm=T),pch=4, col=4, type="b")
lines(x=1:15,y=apply(Haus3,2,mean,na.rm=T),pch=5, col=7, type="b")
lines(x=1:15,y=apply(Haus4,2,mean,na.rm=T),pch=6, col=6, type="b")
legend("topright", legend=c("fused-Lasso","relex fused-Lasso I","relex fused-Lasso II", "relex fused-SCAD I", "relex fused-SCAD II"), pch=c(2,1,4,5,6), col=c(3,2,4,7,6), bty="n")

tau
estimate1
estimate2
estimate3
estimate4
estimate5

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
TP5=matrix(NA,nr=100,nc=15)
TP6=matrix(NA,nr=100,nc=15)
TP7=matrix(NA,nr=100,nc=15)
TP8=matrix(NA,nr=100,nc=15)
FP5=matrix(NA,nr=100,nc=15)
FP6=matrix(NA,nr=100,nc=15)
FP7=matrix(NA,nr=100,nc=15)
FP8=matrix(NA,nr=100,nc=15)
Haus5=matrix(NA,nr=100,nc=15)
Haus6=matrix(NA,nr=100,nc=15)
Haus7=matrix(NA,nr=100,nc=15)
Haus8=matrix(NA,nr=100,nc=15)
set.seed(34)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=rnorm(n,mean=0,sd=0.2)      
  y=signal+e  
  ######relax fused-Lasso I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP5[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP5[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus5[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######relax fused-Lasso II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    TP6[i,length(estimate3)]=TruePositive(estimate3, tau)
    FP6[i,length(estimate3)]=FalsePositive(estimate3, tau)
    Haus6[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
  ######relax fused-SCAD I
  model4=grpreg(X=X1, y=y1, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model4$beta)[2]){
    estimate4=as.vector((which(model4$beta[2:n,j]!=0)))
    if(length(estimate4)<16){
      TP7[i,length(estimate4)]=TruePositive(estimate4, tau)
      FP7[i,length(estimate4)]=FalsePositive(estimate4, tau)
      Haus7[i,length(estimate4)]=hausdorff(estimate4, tau)
    }
  }
  ######relax fused-SCAD II
  model5=grpreg(X=X2, y=y2, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model5$beta)[2]){
    estimate5=as.vector((which(model5$beta[2:n,j]!=0)))
    if(length(estimate5)<16){
      TP8[i,length(estimate5)]=TruePositive(estimate5, tau)
      FP8[i,length(estimate5)]=FalsePositive(estimate5, tau)
      Haus8[i,length(estimate5)]=hausdorff(estimate5, tau)
    }
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
TP9=matrix(NA,nr=100,nc=15)
TP10=matrix(NA,nr=100,nc=15)
TP11=matrix(NA,nr=100,nc=15)
TP12=matrix(NA,nr=100,nc=15)
FP9=matrix(NA,nr=100,nc=15)
FP10=matrix(NA,nr=100,nc=15)
FP11=matrix(NA,nr=100,nc=15)
FP12=matrix(NA,nr=100,nc=15)
Haus9=matrix(NA,nr=100,nc=15)
Haus10=matrix(NA,nr=100,nc=15)
Haus11=matrix(NA,nr=100,nc=15)
Haus12=matrix(NA,nr=100,nc=15)
set.seed(34)
for(i in 1:100){
  tau=((1:J)-1)*n/J+sample(1:(n/J),size=J,replace=TRUE)     #true change points
  beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
  signal=0
  for(j in 1:J) {signal<-signal+beta[j]*(x>tau[j])}
  e=rnorm(n,mean=0,sd=0.2)      
  y=signal+e  
  ######relax fused-Lasso I
  y1=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y
  model2=lars(X1, y=y1, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate2=unlist(model2$actions)[1:j]
    TP9[i,length(estimate2)]=TruePositive(estimate2, tau)
    FP9[i,length(estimate2)]=FalsePositive(estimate2, tau)
    Haus9[i,length(estimate2)]=hausdorff(estimate2, tau)
  }
  ######relax fused-Lasso II
  yy=y-mean(y)
  y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
  model3=lars(X2, y=y2, type="lasso", normalize=FALSE, intercept=FALSE, trace=FALSE, max.steps=20)
  for(j in 1:15){
    estimate3=unlist(model3$actions)[1:j]
    TP10[i,length(estimate3)]=TruePositive(estimate3, tau)
    FP10[i,length(estimate3)]=FalsePositive(estimate3, tau)
    Haus10[i,length(estimate3)]=hausdorff(estimate3, tau)
  }
  ######relax fused-SCAD I
  model4=grpreg(X=X1, y=y1, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model4$beta)[2]){
    estimate4=as.vector((which(model4$beta[2:n,j]!=0)))
    if(length(estimate4)<16){
      TP11[i,length(estimate4)]=TruePositive(estimate4, tau)
      FP11[i,length(estimate4)]=FalsePositive(estimate4, tau)
      Haus11[i,length(estimate4)]=hausdorff(estimate4, tau)
    }
  }
  ######relax fused-SCAD II
  model5=grpreg(X=X2, y=y2, group=1:(n-1), penalty="grSCAD",family="gaussian", gmax=20) 
  for(j in 1:dim(model5$beta)[2]){
    estimate5=as.vector((which(model5$beta[2:n,j]!=0)))
    if(length(estimate5)<16){
      TP12[i,length(estimate5)]=TruePositive(estimate5, tau)
      FP12[i,length(estimate5)]=FalsePositive(estimate5, tau)
      Haus12[i,length(estimate5)]=hausdorff(estimate5, tau)
    }
  }
}

######show Figure 
par(mfrow=c(2,3))
plot(x=1:15,y=apply(TP1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP5,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP9,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP5,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP9,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus1,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus5,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus9,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(TP2,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP6,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP10,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP2,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP6,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP10,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus2,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus6,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus10,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")

par(mfrow=c(2,3))
plot(x=1:15,y=apply(TP3,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP7,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP11,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP3,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP7,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP11,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus3,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus7,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus11,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("h=20","h=40","h=10"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(TP4,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="True positives",main="True positives")
lines(x=1:15,y=apply(TP8,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(TP12,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("bottomright", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(FP4,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,10), xlab = "Number of selected change points",ylab ="False positives",main="False positives")
lines(x=1:15,y=apply(FP8,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(FP12,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topleft", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
plot(x=1:15,y=apply(Haus4,2,mean,na.rm=T),pch=2,col=3, type="b", ylim=c(0,800), xlab = "Number of selected change points",ylab ="Hausdorff distance",main="Hausdorff distance")
lines(x=1:15,y=apply(Haus8,2,mean,na.rm=T), col=2, type="b")
lines(x=1:15,y=apply(Haus12,2,mean,na.rm=T),pch=4, col=4, type="b")
legend("topright", legend=c("kappa=0.5","kappa=1","kappa=0.25"), pch=c(2,1,4), col=c(3,2,4), bty="n")
