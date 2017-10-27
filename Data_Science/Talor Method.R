# The F-function on p.9 of my lecture notes
Ffun <- function(a,b,mu,si,x){
  xs <- (x-mu)/si
  as <- (a-mu)/si
  bs <- (b-mu)/si
  
  f <- function(y) (1/y-1/y^3+3/y^5-15/y^7);
  
  if (as>4 & bs>4){
    F <- 1-(exp((as^2-bs^2)/2)*f(bs) - exp((as^2-xs^2)/2)*f(xs))/
      ( exp((as^2-bs^2)/2)*f(bs) - f(as) )
  } else if (as< -4 & bs< -4) {
    F <- ( exp((bs^2-xs^2)/2)*f(-xs) - exp((bs^2-as^2)/2)*f(-as))/
      ( f(-bs)-exp((bs^2-as^2)/2)*f(-as) )
  } else {
    denom <- pnorm(bs)-pnorm(as)
    if (denom < 0.00001){
      F <- (xs-as)/(bs-as)
    } else {
      F <- (pnorm(xs)-pnorm(as))/denom
    }
  }
  return(F)
}

numSim <- 1000  # number of simulations
Reject_zrand <- rep(0,numSim)
Reject_zmax <- rep(0,numSim)
Reject_taylor <- rep(0,numSim)

Reject_other <- rep(0,numSim)

for (i in 1:numSim){
  # Generate 5 random variables from N(0,1) independently
  Zs <- rnorm(5)
  #Zs <- rnorm(5) + c(2,0,0,0,0)
  
  # Reject zrand if zrand>1.645
  Reject_zrand[i] <- (Zs[sample(5,1)]>qnorm(0.95))
  
  # Reject zmax if zmax>1.645
  jmax <- which(Zs==max(Zs))
  Reject_zmax[i] <- (Zs[jmax]>qnorm(0.95))
  
  # Taylor method: conditional p-value given that jmax is selected
  # Uncomment the following codes when you do part b
  
  # # Compute Gamma and u
  Gamma <- diag(-1,5,5)
  Gamma[,jmax] = Gamma[,jmax] + 1    
  u <- rep(0,5)
  
  # Compute Vup, Vlo
  v <- rep(0,5)
  v[jmax] <- 1
  rho <- (Gamma%*%v)/(sum(v^2))
  if (sum(rho<0)==0){
    Vup <- Inf
  } else {
    P <- c(NULL)
    for (j in 1:length(rho)){
      if (rho[j] < 0){
        P <- c(P,(u[j]-(Gamma%*%Zs)[j]+rho[j]*v%*%Zs)/rho[j])
      }
    }
    Vup <- min(P)
  }
  if (sum(rho>0)==0){
    Vlo <- -Inf
  } else {
    Q <- c(NULL)
    for (j in 1:length(rho)){
      if (rho[j] > 0){
        Q <- c(Q,(u[j]-(Gamma%*%Zs)[j]+rho[j]*v%*%Zs)/rho[j])
      }
    }
    Vlo <- max(Q)
  }
  Reject_taylor[i] <- (Ffun(Vlo,Vup,0,1,Zs[jmax])>0.95)
  
  # Your method for another test that can control type I error at alpha level
  SP <- c(NULL)
  for (j in 1:length(Zs)){
    SP <- c(SP,(1-pnorm(Zs[j])))
  }
  SP <- sort(SP)
  for (j in 1:length(SP)){
    if(SP[j]<(0.05/(5-j+1))){
      Reject_other[i] <- 1
    }
  }
}
print(c(mean(Reject_zrand),mean(Reject_zmax),mean(Reject_taylor),mean(Reject_other
                                                                      
)))