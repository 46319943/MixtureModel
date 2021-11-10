#### Example of an MCMC algorithm for fitting a location mixture of 2 Gaussian components
#### The algorithm is tested using simulated data

## Clear the environment and load required libraries
rm(list=ls())
library(MCMCpack)
set.seed(81196)  # So that results are reproducible

library(readr)
fuses <- read_csv("fuses.csv", col_names = FALSE)
x = fuses$X1
logx = log(x)

KK = 2
n = length(x)

## Initialize the parameters
w = 1/2 # 正常的多于不正常的，所以这个权重可以设小可以点，比如在0.1左右
lambda = n / sum(x) # 不正常的指数分布的参数，可以将n缩小来减少先验概率的有效样本数量。比如设�?0.1*n
mu = mean(logx) #Random cluster centers randomly spread over the support of the data
sigma = sd(logx) #Initial standard deviation

# Plot the initial guess for the density
xx = seq(0, 10, length=100)
yy = w*dexp(xx, lambda) + (1-w)*dlnorm(xx, mu, sigma)
plot(xx, yy, type="l", ylim=c(0, max(yy)), xlab="x", ylab="Initial density")

## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w
# Gamma distribution prior for exponential distribution
alpha = 1
beta = 1
# Normal distribution prior for mean of Log Normal distribution
eta = 0 # Mean 0 for the prior on mu_k
tau = 1 # Standard deviation 5 on the prior for mu_l
# Inverse Gamma distribution prior for varience of Log Normal distribution
dd  = 2
qq  = 1

# Number of iterations of the sampler
rrr   = 6000
burn  = 1000


# Storing the samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
lambda.out = rep(0, rrr)
mu.out    = rep(0, rrr)
sigma.out = rep(0, rrr)
logpost   = rep(0, rrr)

# MCMC iterations
for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = log(w) + dexp(x[i], lambda, log=TRUE)  #Compute the log of the weights
    v[2] = log(1-w) + dlnorm(x[i], mu, sigma, log=TRUE)  #Compute the log of the weights
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))
  
  # Sample the lambda
  lambda = rgamma(1, alpha + sum(cc==1), beta + sum(x[cc==1]))
  
  # Sample the means
  nk = sum(cc==2)
  xsumk = log(prod(x[cc==2]))
  
  tau2.hat = 1/(nk/sigma^2 + 1/tau^2)
  mu.hat = tau2.hat*(xsumk/sigma^2 + eta/tau^2)
  
  mu = rnorm(1, mu.hat, sqrt(tau2.hat))
  
  # Sample the variances
  dd.star = dd + nk/2
  qq.sum = 0
  for(i in 1:n){
    if(cc[i] == 2){
      qq.sum = qq.sum + (log(x[i]) - mu)^2
    }
  }
  qq.star = qq + qq.sum/2
  sigma = sqrt(rinvgamma(1, dd.star, qq.star))
  
  # Store samples
  cc.out[s,] = cc
  w.out[s] = w
  lambda.out[s] = lambda
  mu.out[s] = mu
  sigma.out[s] = sigma
  
  # posterior probability
  
  # P(X | C, theta) * P(C | W)
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w) + dexp(x[i], lambda, log=TRUE)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dlnorm(x[i], mu, sigma, log=TRUE)
    }
  }
  
  # P(W)
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2],log = T)
  
  # P(lambda)
  logpost[s] = logpost[s] + dgamma(lambda, alpha, beta, log = T)
  
  # P(mu)
  logpost[s] = logpost[s] + dnorm(mu, eta, tau, log=T)
  
  # P(sigma)
  logpost[s] = logpost[s] + log(dinvgamma(sigma^2, dd, 1/qq))
  
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

xx = seq(0,10, length=200)
density.posterior = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr-burn)){
  density.posterior[s,] = density.posterior[s,] + w.out[s+burn]*dexp(xx, lambda.out[s+burn]) +
    (1-w.out[s+burn])*dlnorm(xx, mu.out[s+burn], sigma.out[s+burn])
}
density.posterior.m = apply(density.posterior, 2, mean)
density.posterior.lq = apply(density.posterior, 2, quantile, 0.025)
density.posterior.uq = apply(density.posterior, 2, quantile, 0.975)
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(xx, density.posterior.m, type="n",ylim=c(0,max(density.posterior.uq)), xlab="x", ylab="Density")
polygon(c(xx,rev(xx)), c(density.posterior.lq, rev(density.posterior.uq)), col="grey", border="grey")
lines(xx, density.posterior.m, lwd=2)