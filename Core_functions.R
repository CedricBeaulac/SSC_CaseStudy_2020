###########################################
# Important packages
###########################################
require(dlm)
library(forecast)
require(plyr)
require(ggplot2)
require(sinaplot)
require(ggforce)
#require(plot3D)


###########################################
# Fix important variables
###########################################

L= length(total_MW)

# N <- 24*7*4 #Number of training points
# 
# pred <- 24*7*1 #Number of predicted points


###########################################
# Define our own functions
###########################################

bdiag = function(a,b){rbind(cbind(a, matrix(0, nrow=nrow(a), ncol=ncol(b))),
                            cbind(matrix(0, nrow=nrow(b), ncol=ncol(a)), b))}

cr = function(x){(x-mean(x))/sd(x)}

# GetSeries : Extract series at a random time stamp
# inputs : N - Number of training points, pred -Number of predicted points
# output : traning and test points (in a list)
GetSeries <- function(N,pred){
 
  pick <- sample(1:(L-N-pred),1)
  
  y1 <- total_MW[(pick+1):(pick+N)] 
  y2 <- temperature[(pick+1):(pick+N)] 
  
  y1f <- total_MW[(pick+1+N):(pick+N+pred)] 
  y2f <- temperature[(pick+1+N):(pick+N+pred)] 
  
  # y1 <- cr(y1)
  # y2 <- cr(y2)
  # 
  # y1f <- cr(y1f)
  # y2f <- cr(y2f)
  
  return(list(y1,y2,y1f,y2f,pick))
}


# ErrorFunction : MSE
# inputs : y- observations , ypred - predictions
# output : error in R^+
ErrorFunction <- function(y,ypred){
  
  n <- length(y)
  diff <- y-ypred
  squared <- diff^2
  error <- sum(squared)
  
  return(error)
}




##### CV for model 1
ValidationC1 <- function(rs,nt,N,pred){
  

  
  nr <- length(rs)
  
  ErrMatrix <- matrix(data=rep(0,nr),nrow=1)
  
  errors <- rep(0,nt)
  

    for (R in 1:nr){
      
      
      for (i in 1:nt){
        
        series <- GetSeries(N,pred)
        
        prediction <- method1(series,rs[R]) 
        
        errors[i] <- ErrorFunction(prediction$real, prediction$forecast)
        
      }
      
      ErrMatrix[1,R] <- mean(errors)
    }

  
  
  
  return(ErrMatrix)
}


##### CV for model 2
ValidationC2 <- function(qs,rs,nt,N,pred){
  
  nq <- length(qs)
  
  nr <- length(rs)
  
  ErrMatrix <- matrix(data=rep(0,nq*nr),nrow=nq)
  
  errors <- rep(0,nt)
  
  for (Q in 1:nq) {
    for (R in 1:nr){
      
      
      for (i in 1:nt){
        
        series <- GetSeries(N,pred)
        
        prediction <- method2(series,rs[R],qs[Q]) 
        
        errors[i] <- ErrorFunction(prediction$real, prediction$forecast)
        
      }
      
      ErrMatrix[Q,R] <- mean(errors)
    }
  }
  
  
  
  return(ErrMatrix)
}

##### To make CV faster
cv1m3 <- function(N,pred,r,q){
  
  series <- GetSeries(N,pred)
  
  prediction <- method3(series,r,q) 
  
  error <- ErrorFunction(prediction$real, prediction$forecast)
  
  return(error)
  
}


##### CV for model 3
ValidationC3 <- function(qs,rs,nt,N,pred){
  
  nq <- length(qs)
  
  nr <- length(rs)
  
  ErrMatrix <- matrix(data=rep(0,nq*nr),nrow=nq)
  
  errors <- rep(0,nt)
  
  for (Q in 1:nq) {
    for (R in 1:nr){
      
      
      errors <- replicate(nt,cv1m3(N,pred,rs[R],qs[Q]))
      
      ErrMatrix[Q,R] <- mean(errors)
      
      print("One Combination done")
    }
  }
  
  
  
  return(ErrMatrix)
}


##### To make CV faster
cv1m4 <- function(N,pred,r,q){
  
  series <- GetSeries(N,pred)

  prediction <- method4(series,r,q) 

  error <- ErrorFunction(prediction$real, prediction$forecast)
  
  return(error)

}

##### CV for model 4
ValidationC4 <- function(qs,rs,nt,N,pred){
  
  nq <- length(qs)
  
  nr <- length(rs)
  
  ErrMatrix <- matrix(data=rep(0,nq*nr),nrow=nq)
  
  errors <- rep(0,nt)
  
  for (Q in 1:nq) {
    for (R in 1:nr){
      
      errors <- replicate(nt,cv1m4(N,pred,rs[R],qs[Q]))
      
      ErrMatrix[Q,R] <- mean(errors)
    }
  }
  
  
  
  return(ErrMatrix)
}


