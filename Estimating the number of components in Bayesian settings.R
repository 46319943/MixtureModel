EK = 0
alpha = 2
n = 100
for (i in 1:n) {
  EK = EK + alpha / (alpha + i - 1)
}
EK


alpha * log((n + alpha - 1) / alpha) = 2


n = 200
EK = 2
library("nleqslv")
fn <- function(x) {
  return (x * log((n + x - 1) / x) - 2)
}
nleqslv(1, fn)