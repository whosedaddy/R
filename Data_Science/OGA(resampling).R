setwd("C:/Users/User/r")

require("MASS")
require("data.table")
source("OGA_helper.R")

Setting = 1  # set this to 1 for part a, 2 for part b

if (Setting==1){
  # OGA will select 6 predictors out of 6 predictors,
  # which is equivalent to no selection
  # Use this setting for part (a)
  n <- 30; p <- 6; m = 6
}
if (Setting==2){
  # OGA will select 6 predictors out of 15 predictors,
  # Use this setting for part (b)
  n <- 40; p <- 20; m = 6
}

p0 <- 5       # number of non-zero coefficients
alpha <- 0.2  # significance level 
v <- c(rep(1,p0),rep(0,p-p0))

beta <- matrix(rep(0,p),ncol=1)
beta[1:p0] <- rep(1,p0)
theta <- sum(v*beta)
numSimulations <- 500;
Accepts <- rep(0,numSimulations)
AcceptsHB <- rep(0,numSimulations)
AcceptsBoot <- rep(0,numSimulations)
for (i in 1:numSimulations) {
  # Generate X and y for linear regresssion
  X <- matrix(rnorm(n*p),n,p)
  y <- X%*%beta+rnorm(n)
  
  # Do model selection by OGA
  model.select <- OGA(X,y,m)  
  index.select <- model.select$set   # selected index
  betahat <- model.select$beta       # estimated beta with length p
  
  XS <- X[,index.select]
  vs <- v[index.select]
  
  # Method 1: Exact method: assuming normal
  thetahat <- sum(v*betahat)
  sighat <- sqrt(sum((y - X%*%betahat)^2)/(n-m))
  sd <- sighat*sqrt(sum(vs*solve(t(XS)%*%XS,vs)))
  lower <- thetahat - qt(1-alpha/2,df=n-m)*sd
  upper <- thetahat - qt(alpha/2,df=n-m)*sd
  Accepts[i] = (theta>=lower & theta<=upper)
  
  # Method 2: Bootstrap
  Xbeta <- X%*%(solve(t(X)%*%X,t(X)%*%y))
  epshat <- y - Xbeta
  numBoot <- 200
  TS.dist <- rep(0,numBoot)
  for (b in 1:numBoot){
    # Your codes here 
    epshat_b <- sample(epshat,n,T)
    y_b <- Xbeta + epshat_b
	model.select <- OGA(X,y_b,m)
	index.select <- model.select$set
	XS <- X[,index.select]
	vs <- v[index.select]
	betahat_bs <- solve(t(XS)%*%XS,t(XS)%*%y_b)
	thetahat_s <- sum(vs * betahat_bs)
	sig_s <- sqrt(sum((y_b - XS%*%betahat_bs)^2)/(n-m))
	TS.dist[b] <- computeTS(XS,vs,sig_s,thetahat_s,thetahat)
  }
  lower <- thetahat - quantile(TS.dist,1-alpha/2)*sd
  upper <- thetahat - quantile(TS.dist,alpha/2)*sd
  AcceptsBoot[i] = (theta>=lower & theta<=upper)

  # Hybrid resampling
  beta_theta <- ComputeBetaSample(X,y,v,theta)  # assume p<n
  Xbeta_theta <- X%*%beta_theta
  TS.dist <- rep(0,numBoot)
  for (b in 1:numBoot){
    # Your codes here
	epshat_h <- sample(epshat,n,T)
	y_h <- Xbeta_theta + epshat_h
	model.select <- OGA(X,y_h,m)
	index.select <- model.select$set
	XS <- X[,index.select]
	vs <- v[index.select]
	betahat_hs <- solve(t(XS)%*%XS,t(XS)%*%y_h)
	thetahat_h <- sum(vs * betahat_hs)
	sig_h <- sqrt(sum((y_h - XS%*%betahat_hs)^2)/(n-m))
    TS.dist[b] <- computeTS(XS,vs,sig_h,thetahat_h,theta)
  }
  lower <- thetahat - quantile(TS.dist,1-alpha/2)*sd
  upper <- thetahat - quantile(TS.dist,alpha/2)*sd
  
  AcceptsHB[i] = (theta>=lower & theta<=upper)
}
print(c(mean(Accepts),mean(AcceptsBoot),mean(AcceptsHB)))