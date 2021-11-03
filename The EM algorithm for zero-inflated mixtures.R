#### The EM algorithm for zero-inflated mixtures

## Clear the environment and load required libraries
rm(list=ls())
set.seed(81196)    # So that results are reproducible (same simulated data every time)



## read csv
library(readr)
nestsize <- read_csv("nestsize.csv", col_names = FALSE)
x = nestsize$X1

KK = 2
n = length(x)

w = 1/2
lambda = mean(x)

# Plot the initial guess for the density
xx = seq(0,10)
yy = w*c(1,rep(0,length(xx)-1)) + (1-w)*dpois(xx, lambda)
par(mar=c(4,4,2,2)+0.1)
barplot(yy, names.arg=xx, las=1, xlab = "x", ylab="Initial density", 
        border=NA, main="zero-inflated poission mixtures")

s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

zero_mass <- function(x) {
  if(x != 0){
    return(0)
  }
  else{
    return(1)
  }
}

##Checking convergence of the algorithm
while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  # v(i,1)
  for(i in 1:n){
    if(x[i] == 0){
      v[i, 1] = w / (w + (1-w)*dpois(x[i], lambda))
    }
    else{
      v[i, 1] = 0 / (0 + (1-w)*dpois(x[i], lambda))
    }
  }
  # v(i,2)
  for(i in 1:n){
    if(x[i] == 0){
      v[i, 2] = (1-w)*dpois(x[i], lambda) / (w + (1-w)*dpois(x[i], lambda))
    }
    else{
      v[i, 2] = (1-w)*dpois(x[i], lambda) / (0 + (1-w)*dpois(x[i], lambda))
    }
  }
  
  ## M step
  # Weights
  w = mean(v[,1])
  mu = rep(0, KK)
  for(k in 1:KK){
    for(i in 1:n){
      mu[k]    = mu[k] + v[i,k]*x[i]
    }
    mu[k] = mu[k]/sum(v[,k])
  }
  
  # Standard deviations
  sigma = 0
  for(i in 1:n){
    for(k in 1:KK){
      sigma = sigma + v[i,k]*(x[i] - mu[k])^2
    }
  }
  sigma = sqrt(sigma/sum(v))
  
  ##Check convergence
  QQn = 0
  for(i in 1:n){
    QQn = QQn + v[i,1]*(log(w) + dnorm(x[i], mu[1], sigma, log=TRUE)) +
      v[i,2]*(log(1-w) + dnorm(x[i], mu[2], sigma, log=TRUE))
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  
  #Plot current estimate over data
  layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
  par(mar=c(3.1,4.1,0.5,0.5))
  plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)
  
  par(mar=c(5,4,1.5,0.5))
  xx = seq(-8,11,length=200)
  yy = w*dnorm(xx, mu[1], sigma) + (1-w)*dnorm(xx, mu[2], sigma)
  plot(xx, yy, type="l", ylim=c(0, max(c(yy,yy.true))), main=paste("s =",s,"   Q =", round(QQ.out[s],4)), lwd=2, col="red", lty=2, xlab="x", ylab="Density")
  lines(xx.true, yy.true, lwd=2)
  points(x, rep(0,n), col=cc)
  legend(6,0.22,c("Truth","Estimate"),col=c("black","red"), lty=c(1,2))
}


#Plot final estimate over data
layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
par(mar=c(3.1,4.1,0.5,0.5))
plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)

par(mar=c(5,4,1.5,0.5))
xx = seq(-8,11,length=200)
yy = w*dnorm(xx, mu[1], sigma) + (1-w)*dnorm(xx, mu[2], sigma)
plot(xx, yy, type="l", ylim=c(0, max(c(yy,yy.true))), main=paste("s =",s,"   Q =", round(QQ.out[s],4)), lwd=2, col="red", lty=2, xlab="x", ylab="Density")
lines(xx.true, yy.true, lwd=2)
points(x, rep(0,n), col=cc)
legend(6,0.22,c("Truth","Estimate"),col=c("black","red"), lty=c(1,2), bty="n")


