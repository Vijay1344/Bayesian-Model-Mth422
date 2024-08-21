Y      <- c(64,13,33,18,30,20)

post   <- function(Y,t,beta,pri.sd=10){
  mn <- exp(beta[1] + t*beta[2])
  l  <- prod(dpois(Y,mn))
  p  <- prod(dnorm(beta,0,pri.sd))
  return(l*p)}

MCMC <- function(Y, beta.init, iters,cansd){
  
  # chain initiation
  beta <- beta.init
  # define chains
  beta.chain <- matrix(NA, iters,2)
  
  t <- 1:length(Y)
  # start MCMC
  for(iter in 1:iters){
    for(j in 1:2){
      can    <- beta
      can[j] <- rnorm(1,beta[j],cansd[j])
      R      <- post(Y,t,can)/post(Y,t,beta)
      if(runif(1)<R){
        beta <- can
      }
    }
    beta.chain[iter,] <- beta
  }
  
  # return chains
  out <- list(beta.chain = beta.chain)
  return(out)}  


MCMC.out <- MCMC(Y, beta.init = c(2,0),cansd = c(0.2,0.05), iters = 30000)
plot(MCMC.out$beta.chain[,1],type = "l",ylab = expression(alpha), xlab = "Iteration")
plot(MCMC.out$beta.chain[,2],type = "l",ylab = expression(beta), xlab = "Iteration")


### acceptance ratio
S <- 30000
acc_rate <- colMeans(MCMC.out$beta.chain[-1,]!=MCMC.out$beta.chain[-S,])
acc_rate
