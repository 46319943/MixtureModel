# Generate n observations from a mixture of three Poisson 
# distributions
n     = 200           # Size of the sample to be generated
w     = c(0.7, 0.2, 0.1)  # Weights
lambda = c(1, 2, 6)
cc    = sample(1:3, n, replace=T, prob=w)
x     = rpois(n, lambda[cc])

# Plot f(x) along with the observations 
# just sampled
xx = seq(0, 12)
yy = w[1]*dpois(xx, lambda[1]) + w[2]*dpois(xx, lambda[2]) + w[3]*dpois(xx, lambda[3])
par(mar=c(4,4,1,1)+0.1)
barplot(yy)
empfreq = table(factor(x, levels=seq(0, max(x))))/n
barplot(empfreq)
hist(x, right=F)
