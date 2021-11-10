# 搞不懂这个分布。
# 最后还是就当它是一个普通的分布求的后验概率

mean = 3
sd = 1

x = 10
lnx = log(x)

dlnorm(x, mean, sd)
dnorm(lnx, mean, sd)