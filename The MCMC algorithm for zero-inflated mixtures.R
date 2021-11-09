#### Example of an MCMC algorithm for fitting a location mixture of 2 Gaussian components
#### The algorithm is tested using simulated data

## Clear the environment and load required libraries
rm(list=ls())
library(MCMCpack)
set.seed(81196)  # So that results are reproducible


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

## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w

# Gamma for Poisson prior
alpha = 1 
beta = 1

# Number of iterations of the sampler
rrr   = 6000
burn  = 1000


# Storing the samples
cc.out = array(0, dim=c(rrr, n))
w.out = rep(0, rrr)
lambda.out = array(0, rrr)
logpost = rep(0, rrr)

# MCMC iterations
for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    
    if(x[i] == 0){
      v[1] = w
      v[2] = (1-w)*dpois(x[i], lambda)
    }
    else{
      v[1] = 0
      v[2] = 1
    }
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))
  
  # Sample the lambda
  lambda = rgamma(1, alpha + sum(x[cc==2]), beta + sum(cc==2))
  
  # Store samples
  cc.out[s,] = cc
  w.out[s] = w
  lambda.out[s] = lambda
  
  # posterior probability
  
  # P(X | C, theta) * P(C | W)
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dpois(x[i], lambda, log=TRUE)
    }
  }
  
  # P(w)
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2],log = T)
  
  # P(theta)
  logpost[s] = logpost[s] + dgamma(lambda, 1, 1, log = T)
  
  
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

xx = seq(0,10)
density.posterior = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr-burn)){
  density.posterior[s,] = density.posterior[s,] + 
    w.out[s+burn]*c(1,rep(0,length(xx)-1)) + 
    (1-w.out[s+burn])*dpois(xx,lambda.out[s+burn])
}
density.posterior.m = apply(density.posterior , 2, mean)
density.posterior.lq = apply(density.posterior, 2, quantile, 0.025)
density.posterior.uq = apply(density.posterior, 2, quantile, 0.975)
par(mfrow=c(1,1))

par(mar=c(4,4,2,2)+0.1)
barplot(density.posterior.m, names.arg=xx, las=1, xlab = "x", ylab="Initial density", 
        border=NA, main="zero-inflated poission mixtures")

polygon(c(xx,rev(xx)), c(density.posterior.lq, rev(density.posterior.uq)), col="grey", border="black")
lines(xx, density.posterior.m, lwd=2)
