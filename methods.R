#### Method 1 : Univariate LSSM model
#### LSSM stochastic: local level (mean), deterministic: cyc day + cyc week + cyc year
#### Parameters: r signal to noise ratio sigma_mu/sigma_y

#* is a line that is useless to repeat in a CV run; should be taken out of loop. Check this later
## NB. should be negligible vs. time to fit model


method1 <- function(series, r){
  
  
  ### useful variables
  y_1toN <- series[[1]]
  pred <- length(series[[3]]) #*
  ########
    
    
  ### building model
  nbHarm = 3 #*
  dlmModel <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*7,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) #*
  dlmModel$V <- 1 #* 
  dlmModel$W[1,1] <- r
  ##############                    
  
  
  yFilt <- dlmFilter(y_1toN,dlmModel)
  yFore <- dlmForecast(yFilt,pred)
  
  
  return(list(real = series[[3]], forecast = as.numeric(yFore$f)))
}


#### Method 3 : Our model presented in 2020 SSC Case study (multivariate)
#### bivariate LSSM stochastic: 2 correlated local level (mean), deterministic: 2 cyc day + cyc week + 2 cyc year
#### q is negative in the winter and positive in the summer
#### Parameters: r signal to noise ratio sigma_mu/sigma_y, q>0 level correlation


method3 <- function(series, r, q){
  
  
  ### useful variables
  N <- length(series[[1]]) #*
  pred <- length(series[[3]]) #* 
  y_1toN <- cr(series[[1]])
  x <- cr(c(series[[2]],series[[4]]))
  x_1toN <- x[1:N]
  xFore <- x[(N+1):(N+pred)]
  pick <- series[[5]]
  year <- c(rep(-q,151*24),rep(q,122*24),rep(-q,92*24)) 
  X <- rep(year,20) 
  X <- X[(pick+1):length(X)] 
  ########
  
  #* 
  nbHarm = 3 
  dlmModel1 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*7,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  dlmModel2 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  
  dlmModel <- dlmModel1
  
  dlmModel$FF <- bdiag(FF(dlmModel1),FF(dlmModel2))
  
  dlmModel$V <- diag(2)
  
  dlmModel$GG <- bdiag(GG(dlmModel1),GG(dlmModel2))
  
  dlmModel$W <- bdiag(W(dlmModel1),W(dlmModel2))
  
  dlmModel$m0 <- c(m0(dlmModel1),m0(dlmModel2))
  dlmModel$C0 <- bdiag(C0(dlmModel1),C0(dlmModel2))
  
  secMean <- 2 + nbHarm*2*3
  #nbState <- length(dlmModel$m0)
  #*
  
  
  
  dlmModel$W[1,1] <- r
  dlmModel$W[secMean,secMean] <- r
  
  dlmModel$JW <- dlmModel$W
  
  dlmModel$JW[1,1] <- 0
  dlmModel$JW[secMean,secMean] <- 0
  
  dlmModel$JW[1,secMean] <- 1
  dlmModel$JW[secMean,1] <- 1
  
  dlmModel$X <- as.matrix(X)
  
  ptm <- proc.time()
  
  
  yFilt <- dlmFilter(cbind(c(y_1toN,rep(NA,pred)),x), dlmModel)
  proc.time() - ptm
  
  yFore <- as.numeric(t(as.matrix(dlmModel$FF[1,]))%*%t(yFilt$m[(N+2):(N+pred+1),]))*sd(series[[1]]) + mean(series[[1]])
  
  
  return(list(real = series[[3]], forecast = yFore))
}



###########################################
# Alternative models not presented
###########################################


#### Method 2
#### bivariate LSSM stochastic: 2 correlated local level (mean), deterministic: 2 cyc day + cyc week + 2 cyc year
#### Parameters: r signal to noise ratio sigma_mu/sigma_y, q level correlation


method2 <- function(series, r, q){
  
  
  ### useful variables
  N <- length(series[[1]]) #*
  pred <- length(series[[3]]) #* 
  y_1toN <- cr(series[[1]])
  x <- cr(c(series[[2]],series[[4]]))
  x_1toN <- x[1:N]
  xFore <- x[(N+1):(N+pred)]
  
  ########
  
  #* 
  nbHarm = 3 
  dlmModel1 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*7,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  dlmModel2 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  
  dlmModel <- dlmModel1
  
  dlmModel$FF <- bdiag(FF(dlmModel1),FF(dlmModel2))
  
  dlmModel$V <- diag(2)
  
  dlmModel$GG <- bdiag(GG(dlmModel1),GG(dlmModel2))
  
  dlmModel$W <- bdiag(W(dlmModel1),W(dlmModel2))
  
  dlmModel$m0 <- c(m0(dlmModel1),m0(dlmModel2))
  dlmModel$C0 <- bdiag(C0(dlmModel1),C0(dlmModel2))
  
  secMean <- 2 + nbHarm*2*3
  #nbState <- length(dlmModel$m0)
  #*
  
  dlmModel$W[1,1] <- r
  dlmModel$W[secMean,secMean] <- r
  
  
  dlmModel$W[1,secMean] <- q*r
  dlmModel$W[secMean,1] <- q*r
  
  yFilt <- dlmFilter(cbind(c(y_1toN,rep(NA,pred)),x), dlmModel)
  
  
  yFore <- as.numeric(t(as.matrix(dlmModel$FF[1,]))%*%t(yFilt$m[(N+2):(N+pred+1),]))*sd(series[[1]]) + mean(series[[1]])
  
  
  return(list(real = series[[3]], forecast = yFore))
}



#### Method 4 
#### bivariate LSSM stochastic: 2 correlated local level (mean), deterministic: 2 cyc day + cyc week + 2 cyc year
#### temprature is encoded as |temperature-21|
#### Parameters: r signal to noise ratio sigma_mu/sigma_y, q level correlation


method4 <- function(series, r, q){
  
  
  ### useful variables
  N <- length(series[[1]]) #*
  pred <- length(series[[3]]) #* 
  y_1toN <- cr(series[[1]])
  x <- cr(abs(c(series[[2]],series[[4]])-21))
  x_1toN <- x[1:N]
  xFore <- x[(N+1):(N+pred)]
  
  ########
  
  #* 
  nbHarm = 3 
  dlmModel1 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*7,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  dlmModel2 <- dlmModPoly(order=1)  + dlmModTrig(s=24,q=nbHarm) + dlmModTrig(s=24*365,q=nbHarm) 
  
  dlmModel <- dlmModel1
  
  dlmModel$FF <- bdiag(FF(dlmModel1),FF(dlmModel2))
  
  dlmModel$V <- diag(2)
  
  dlmModel$GG <- bdiag(GG(dlmModel1),GG(dlmModel2))
  
  dlmModel$W <- bdiag(W(dlmModel1),W(dlmModel2))
  
  dlmModel$m0 <- c(m0(dlmModel1),m0(dlmModel2))
  dlmModel$C0 <- bdiag(C0(dlmModel1),C0(dlmModel2))
  
  secMean <- 2 + nbHarm*2*3
  #nbState <- length(dlmModel$m0)
  #*
  
  dlmModel$W[1,1] <- r
  dlmModel$W[secMean,secMean] <- r
  
  
  dlmModel$W[1,secMean] <- q*r
  dlmModel$W[secMean,1] <- q*r
  
  yFilt <- dlmFilter(cbind(c(y_1toN,rep(NA,pred)),x), dlmModel)
  
  
  yFore <- as.numeric(t(as.matrix(dlmModel$FF[1,]))%*%t(yFilt$m[(N+2):(N+pred+1),]))*sd(series[[1]]) + mean(series[[1]])
  
  
  return(list(real = series[[3]], forecast = yFore))
}






