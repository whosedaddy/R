# You should not change any codes in this file

# This function is to help you compute T* in 
# Bootstrap and Hybrid resampling 
# You may write your own T* if you disagree or don't understand what 
# this function is doing
computeTS <- function(X,v,sig,thetahat,theta){
  denom <- sig*sqrt(t(v)%*%solve(t(X)%*%X,v))
  return((thetahat-theta)/denom)
}

# This function is to compute betahat under H0 for 
# Hybrid resampling
# The output beta will be of length p
ComputeBetaSample <- function(X,y,v,theta){
  KKT.matrix <- rbind( cbind(t(X)%*%X, v), c(t(v),0))
  KKT.vec <- c(t(X)%*%y, theta)
  beta_theta <- solve(KKT.matrix, KKT.vec)
  return(beta_theta[1:(length(beta_theta)-1)])
}

# This OGA function selects m predictors out of p predictors
# No need to understand how OGA do the selection
# Just need to know that it generate a selected set 
# and a beta estimate
OGA <- function(X,y,m){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  U <- y
  j <- rep(0,m)
  Xnorms <- sqrt(apply(X^2,2,sum))
  Qk <- matrix(rep(0,n*m),n,m)
  Rk <- matrix(rep(0,m*m),m,m)
  beta <- rep(0,p)
  BetaTemp <- rep(0,m)
  sigmahat2 <- rep(0,m+1) #include no selection
  
  for (k in 1:m){
    sigmahat2[k] = sum(U^2)
    
    C <- abs((t(U)%*%X))/Xnorms
    if (k>1){
      C[j[1:k-1]] <- 0
    }
    j[k] <- which(C==max(C))[1]
    
    Xk = X[,j[k]]
    if (k==1){
      rk <- sum(Xk^2)
      qk <- Xk/rk
      Qk[,1] <- qk
      Rk[1,1] <- rk
    } else {
      w <- t(Qk[,(1:k-1)])%*%Xk
      rq <- Xk - Qk[,(1:k-1)]%*%w
      rk <- sum(rq^2)
      qk <- rq/rk
      Qk[,k] <- qk
      Rk[k,k] <- rk
      Rk[1:(k-1),k] <- w
    }
    
    BetaTemp[k] = sum(qk*U)
    U = U - qk*BetaTemp[k] 
  }
  sigmahat2[m+1] <- sum(U^2)
  
  XS <- X[,j]
  beta[j] <- solve(t(XS)%*%XS,t(XS)%*%y)
  
  return(list(beta=beta,set=j))
}