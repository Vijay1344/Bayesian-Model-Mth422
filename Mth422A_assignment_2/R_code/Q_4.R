##########4
#Generate simulated data
set.seed(123) 
T <- 10
rho <- 0.5
X <- numeric(T)
X[1] <- rnorm(1, 0, 1)
for (t in 2:T) {
  X[t] <- rnorm(1, mean = rho * X[t-1], sd = sqrt(1 - rho^2))
}

density.rho.given.data <- function(rho)
{
  X.t <- 0
  for (t in 2:T)
  {
    X.t <- (X.t +  sum((X[t])- rho*(X[t-1])))
  }
  e.fun.val <- exp(-X.t/(2*(1-rho^2)))
  
  densty <- (2*pi)^(-T/2) *(1-rho^2)^( (1-T)/2) * exp( (-X[1]^2)/2 ) *e.fun.val
  return(densty)
}
p <- seq(-1,1, 0.005)
plot(p, density.rho.given.data(p), type = "l")



#### max val occur at rho = 0 by viewing fig
c <- density.rho.given.data(0)

accept.reject <- function()
{
  accept <- 0
  while (accept == 0) {
    pros <- runif(1, min = -1, max = 1) ## proposal density U(-1,1)
    U <- runif(1)
    ratio <- density.rho.given.data(pros)/(c*dunif(pros, min = -1, max = 1))
    if(U <= ratio)
    {
      accept <- 1
      return(pros)
    }
    
  }
}


N <- 1e5
samp <- numeric(length = N)

counts <- numeric(length = N)
for(i in 1:N)
{
  samp[i] <- accept.reject()
}
plot(density(samp), main = " kernel density estimate of rho|X1,..XT", xlab = "rho")
