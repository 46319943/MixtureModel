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
epsilon = 0.00001

par(mfrow=c(1,1))
par(mar=c(2,2,2,2)+0.1)
colscale = c("black","red")
pairs(
  banknote.training[,-1], 
  # 颜色
  col=colscale[banknote.training.labels],
  # 形状
  pch=training_label_numeric
  )


# Initialize the parameters of the algorithm
set.seed(63252)
w   = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu  = rmvnorm(KK, apply(x,2,mean), var(x))   #Cluster centers randomly spread over the support of the data
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/KK  
Sigma[2,,] = var(x)/KK
Sigma[3,,] = var(x)/KK

sw     = FALSE
KL     = -Inf
KL.out = NULL
s      = 0



while(!sw){
  ## E step
  v = array(0, dim=c(n+m,KK))
  for(k in 1:KK){  #Compute the log of the weights
    v[1:n,k] = ifelse(training_label_numeric==k,0,-Inf)  # Training set
    v[(n+1):(n+m),k] = log(w[k]) + mvtnorm::dmvnorm(x[(n+1):(n+m),], mu[k,], Sigma[k,,],log=TRUE)  # Test set
  }
  for(i in 1:(n+m)){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  }
  
  ## M step
  w = apply(v,2,mean)
  mu = matrix(0, nrow=KK, ncol=p)
  for(k in 1:KK){
    for(i in 1:(n+m)){
      mu[k,]    = mu[k,] + v[i,k]*x[i,]
    }
    mu[k,] = mu[k,]/sum(v[,k])
  }
  Sigma = array(0,dim=c(KK,p,p))
  for(k in 1:KK){
    for(i in 1:(n+m)){
      Sigma[k,,] = Sigma[k,,] + v[i,k]*(x[i,] - mu[k,])%*%t(x[i,] - mu[k,])
    }
    Sigma[k,,] = Sigma[k,,]/sum(v[,k])
  }
  
  ##Check convergence
  KLn = 0
  for(i in 1:(n+m)){
    for(k in 1:KK){
      KLn = KLn + v[i,k]*(log(w[k]) + mvtnorm::dmvnorm(x[i,],mu[k,],Sigma[k,,],log=TRUE))
    }
  }
  if(abs(KLn-KL)/abs(KLn)<epsilon){
    sw=TRUE
  }
  KL = KLn
  KL.out = c(KL.out, KL)
  s = s + 1
  print(paste(s, KLn))
}

## Predicted labels
apply(v, 1, which.max)[(n+1):(n+m)]
## True labels
test_label_numeric
## Comparison
apply(v, 1, which.max)[(n+1):(n+m)] == test_label_numeric
sum(!(apply(v, 1, which.max)[(n+1):(n+m)] == test_label_numeric)) # One error

# Using the qda and lda functions in R
# qda
modqda = qda(grouping=training_label_numeric, x=banknote.training, method="mle")
ccpredqda = predict(modqda,newdata=banknote.test)
sum(!(ccpredqda$class == test_label_numeric)) # 3 error

# lda
modlda = lda(grouping=training_label_numeric, x=banknote.training, method="mle")
ccpredlda = predict(modlda,newdata=banknote.test)
sum(!(ccpredlda$class == test_label_numeric)) # 1 error