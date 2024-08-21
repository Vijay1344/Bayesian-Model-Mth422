set.seed(1)
n = 15
m = 10
lambda.1 <- 2
lambda.2 <- 2.5
a <- 0.1
b <- 0.1
X <- rpois(n,lambda = lambda.1)
Y <- rpois(m, lambda = lambda.2)

lambda.1.samples <- rgamma(1e4, shape = sum(X) + a, rate = m+b)
lambda.2.samples <- rgamma(1e4, shape = sum(Y) + a, rate = n + b)
theta.samples <- lambda.1.samples/(lambda.1.samples + lambda.2.samples)
#  95% HPD CI based on Monte Carlo samples
opt.x.s <- optimize(function(x){
  l <- quantile(theta.samples, c(x,x + 0.95))
  return(l[2]-l[1])
}, lower = 0,upper = 0.05)$minimum


posCI95.HPD.s <- quantile(theta.samples, c(opt.x.s,opt.x.s + 0.95))

print("posCI95.HPD.s")
posCI95.HPD.s