# Generate n observations from a mixture of two Gaussian 
# distributions
n     = 100           # Size of the sample to be generated
w     = c(0.3, 0.25, 0.25, 0.2)  # Weights
mu    = c(1, 4, 7, 10)      # Means
cc    = sample(1:4, n, replace=T, prob=w)
x     = rexp(n, 1 / mu[cc])

approximate_mean = mean(x)
approximate_variance = var(x)