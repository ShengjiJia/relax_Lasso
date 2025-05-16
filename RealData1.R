library(seqbbs)   #is available at https://github.com/metalhelix/seqbbs
library(lars)
library(grpreg)
library(KernSmooth)
library(DNAcopy)
library(imputeTS)

#################################define some functions
estimateSigma<-function (Y, h = 10) {  
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

localMax<-function (y, span = 5) {  
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

#################################real data 3
test_filename <- system.file("extdata", "paper.txt", package="seqbbs")
ratios <- read.table(test_filename, header = FALSE)
y=ratios$V1               
n=length(y)     
x=1:n
x1=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) {x1[i,(i+1):n]=0}                   
##############################method1 CBS
set.seed(18)
CBS=DNAcopy::segment(CNA(as.matrix(y), chrom=rep(1,n), maploc=1:n))
estimate1=sort(unique(CBS$output[,4]))
estimate1=estimate1[-length(estimate1)]
##############################method2 SaRa
D=abs(localDiagnostic(y, h=10))
index=localMax(D,span=10)
lambda=2*sqrt(2/10)*mad(diff(y))/sqrt(2)        #use robust estimator of sigma
candidate=index[which(D[index]>lambda)]
candidate1=candidate[order(D[candidate],decreasing=TRUE)]          #ranking
bic=rep(0,length(candidate1))
for(j in 1:length(candidate1)){
  model=lm(y~x1[,candidate1[1:j]])
  bic[j]=n*log(summary(model)$sigma)+j*log(n)
}
j1=which.min(bic)
estimate2=sort(candidate1[1:j1])
##############################method3 fused Lasso  
group=1:n-1          
model1=grpreg(X=x1[,2:n], y=y, group=group[2:n], penalty="grLasso",family="gaussian", gmax=100) 
opt=which.min(log(model1$deviance/n)+model1$df*0.2*log(n)*log(log(n))/n)   #optimal lambda
estimate=as.vector((which(model1$beta[,opt]!=0)[-1]))-1
estimate3=refine(estimate,error=4)
##############################method4 Proposed I
h=20       
xx1=x1
for(i in 2:n) {xx1[,i]=x1[,i]-locpoly(x,x1[,i],kernel="epanech",bandwidth=h,gridsize=n)$y }      
y0=y-locpoly(x,y,kernel="epanech",bandwidth=h,gridsize=n)$y           
model2=grpreg(X=xx1[,2:length(y0)], y=y0, group=group[2:length(y0)], penalty="grLasso",family="gaussian", gmax=100) 
opt=which.min(log(model2$deviance/n)+model2$df*0.2*log(n)*log(log(n))/n)   #optimal lambda
estimate=as.vector((which(model2$beta[,opt]!=0)[-1]))-1
estimate4=refine(estimate,error=4)
##############################method5 Proposed II
lambda=1
X=x1[,2:n]
XX=X                                   #centered X
for(i in 1:(n-1)) {XX[,i]=X[,i]-mean(X[,i])}
R=XX%*%solve(t(XX)%*%XX+n*lambda*diag(1,n-1))%*%t(XX)
eig=eigen(diag(1,n)-R)
X2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%XX       #matrix (I-R)^{1/2}(I-11/n)X
yy=y-mean(y)
y2=eig$vectors%*%sqrt(diag(eig$values))%*%t(eig$vectors)%*%yy
model3=grpreg(X=X2, y=y2, group=group[2:length(y2)], penalty="grLasso",family="gaussian", gmax=100)
opt=which.min(log(model3$deviance/n)+model3$df*0.2*log(n)*log(log(n))/n)   #optimal lambda
estimate=as.vector((which(model3$beta[,opt]!=0)[-1]))-1
estimate5=refine(estimate,error=4)
##############################show figures
par(mfrow=c(2,3))
plot(x, y, xlim=c(300,800), ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="CBS", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate1)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate1,n)[i]+1):(c(0,estimate1,n)[i+1])]), c(0,estimate1,n)[i+1]-c(0,estimate1,n)[i]))
}
residual=y-fitted
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="SaRa", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate2)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate2,n)[i]+1):(c(0,estimate2,n)[i+1])]), c(0,estimate2,n)[i+1]-c(0,estimate2,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="fused-Lasso", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate3)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate3,n)[i]+1):(c(0,estimate3,n)[i+1])]), c(0,estimate3,n)[i+1]-c(0,estimate3,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="Proposed I", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate4)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate4,n)[i]+1):(c(0,estimate4,n)[i+1])]), c(0,estimate4,n)[i+1]-c(0,estimate4,n)[i]))
}
lines(fitted,col=2,lwd=2)
plot(x, y, xlim=c(300,800),  ylim=c(0,2.5), xlab="locations", ylab="Log 2 ratio", main="Proposed II", pch=20, col=8)
fitted=NULL
for(i in 1:(length(estimate5)+1)){
  fitted=c(fitted,rep(mean(y[(c(0,estimate5,n)[i]+1):(c(0,estimate5,n)[i+1])]), c(0,estimate5,n)[i+1]-c(0,estimate5,n)[i]))
}
lines(fitted,col=2,lwd=2)
pacf(residual,lag.max=20,main="PACF")

print(c(length(estimate1),length(estimate2),length(estimate3),length(estimate4),length(estimate5)))
