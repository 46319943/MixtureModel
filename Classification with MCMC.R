# Semi-supervised, quadratic discriminat analysis 
### Loading data and setting up global variables
load("banknoteclassification.Rdata")
training_label_numeric = as.numeric(
  banknote.training.labels
)
test_label_numeric = as.numeric(
  banknote.test.labels
)
library(MASS)
library(mvtnorm)
n = dim(banknote.training)[1]  # Size of the training set
m = dim(banknote.test)[1]      # Size of the test set
x = rbind(as.matrix(banknote.training), as.matrix(banknote.test))   # Create dataset of observations, first n belong to the training set, and the rest belong to the test set
p       = dim(x)[2]              # Number of features
KK      = 2 # 分类数量

# Using MCMC
library(ellipse)
library(MCMCpack)

## Initialize the parameters
set.seed(63252)
w          = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu         = rmvnorm(KK, apply(x,2,mean), var(x))   #RandomCluster centers randomly spread over the support of the data
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/KK  
Sigma[2,,] = var(x)/KK
cc         = sample(1:KK, n, replace=TRUE, prob=w)

par(mfrow=c(1,1))
# 只选取了前两个特征
plot(x[,1], x[,2], 
     col=append(training_label_numeric, test_label_numeric), 
     xlab=expression(x[1]), ylab=expression(x[2]))
for(k in 1:KK){
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.50), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.82), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.95), col="grey", lty=2, lwd=2)
}
title(main="Initial estimate + Observations")

# Priors
aa = rep(1, KK)
dd = apply(x,2,mean)
DD = 10*var(x)
nu = p
SS = var(x)/KK

# Number of iteration of the sampler
rrr = 1000
burn = 100

# Storing the samples
cc.out    = array(0, dim=c(rrr, n+m))
w.out     = array(0, dim=c(rrr, KK))
mu.out    = array(0, dim=c(rrr, KK, p))
Sigma.out = array(0, dim=c(rrr, KK, p, p))
logpost   = rep(0, rrr)

for(s in 1:rrr){
  # Sample the indicators
  for(i in 1:(n+m)){
    v = rep(0,KK)
    for(k in 1:KK){
      # training set
      if(i <= n){
        v[k] = ifelse(training_label_numeric[i]==k, 0, -Inf)
      }
      # test set
      else{
        v[k] = log(w[k]) + mvtnorm::dmvnorm(x[i,], mu[k,], Sigma[k,,], log=TRUE)  #Compute the log of the weights
      }
    }
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = as.vector(rdirichlet(1, aa + tabulate(cc)))
  
  # Sample the means
  DD.st = matrix(0, nrow=p, ncol=p)
  for(k in 1:KK){
    mk    = sum(cc==k)
    xsumk = apply(x[cc==k,], 2, sum)
    DD.st = solve(mk*solve(Sigma[k,,]) + solve(DD))
    dd.st = DD.st%*%(solve(Sigma[k,,])%*%xsumk + solve(DD)%*%dd)
    mu[k,] = as.vector(rmvnorm(1,dd.st,DD.st))
  }
  
  # Sample the variances
  xcensumk = array(0, dim=c(KK,p,p))
  # sample size is n+m, not just n
  for(i in 1:n+m){
    xcensumk[cc[i],,] = xcensumk[cc[i],,] + (x[i,] - mu[cc[i],])%*%t(x[i,] - mu[cc[i],])
  }
  for(k in 1:KK){
    Sigma[k,,] = riwish(nu + sum(cc==k), SS + xcensumk[k,,])
  }
  
  # Store samples
  cc.out[s,]      = cc
  w.out[s,]       = w
  mu.out[s,,]     = mu
  Sigma.out[s,,,] = Sigma
  for(i in 1:n){
    logpost[s] = logpost[s] + log(w[cc[i]]) + mvtnorm::dmvnorm(x[i,], mu[cc[i],], Sigma[cc[i],,], log=TRUE)
  }
  logpost[s] = logpost[s] + ddirichlet(w, aa)
  for(k in 1:KK){
    logpost[s] = logpost[s] + mvtnorm::dmvnorm(mu[k,], dd, DD, log=TRUE)
    logpost[s] = logpost[s] + log(diwish(Sigma[k,,], nu, SS))
  }
  
  if(s/250==floor(s/250)){
    print(paste("s = ", s))
  }
  
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

## Plot the density estimate for the last iteration of the MCMC
par(mfrow=c(1,1))
par(mar=c(4,4,2,1)+0.1)
pair_feature_1 = 3
pare_feature_2 = 6
plot(x[,pair_feature_1], x[,pare_feature_2], 
     col=append(training_label_numeric, test_label_numeric),
     main=paste("s =",s,"   logpost =", round(logpost[s],4)),
     xlab=expression(x[pair_feature_1]), ylab=expression(x[pare_feature_2]))

# TODO: ellipse which参数不懂，怎么指定用哪个维度进行绘制？
for(k in 1:KK){
  lines(
    ellipse(
      x=Sigma[k,c(pair_feature_1,pare_feature_2),c(pair_feature_1,pare_feature_2)], 
      centre=mu[k,c(pair_feature_1,pare_feature_2)], 
      level=0.50, 
    )
    ,col="grey", lty=2, lwd=2
  )
  lines(
    ellipse(
      x=Sigma[k,c(pair_feature_1,pare_feature_2),c(pair_feature_1,pare_feature_2)], 
      centre=mu[k,c(pair_feature_1,pare_feature_2)], 
      level=0.82, 
      )
    ,col="grey", lty=2, lwd=2
    )
  lines(
    ellipse(
      x=Sigma[k,c(pair_feature_1,pare_feature_2),c(pair_feature_1,pare_feature_2)], 
      centre=mu[k,c(pair_feature_1,pare_feature_2)], 
      level=0.95, 
    )
    ,col="grey", lty=2, lwd=2
  )
}

# TODO: 需要使用cc的概率去判断属于哪一类

# drop burn and average the cc.out result, 0 error
sum(
  !round(
    colMeans(
      cc.out[burn:rrr,(n+1):(n+m)]
    )
  )
)

# use max log-likelihood iteration, 3 error
max_log_index = which.max(logpost)
sum(
  !(cc.out[max_log_index,][(n+1):(n+m)] == test_label_numeric)
  )

# use the last iteration, 1 error
sum(
  !(cc.out[rrr,][(n+1):(n+m)] == test_label_numeric)
)

# qda
modqda = qda(grouping=training_label_numeric, x=banknote.training, method="mle")
ccpredqda = predict(modqda,newdata=banknote.test)
sum(!(ccpredqda$class == test_label_numeric)) # 3 error
