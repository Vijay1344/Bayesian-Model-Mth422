library(rjags)
Y <- 10
a <- 1
b <- 1
lambda <- 10

data <- list(Y = Y, a = a,b = b, lambda = lambda)

model_string <- textConnection("model{
    # liklihood
    Y ~ dbin(p,n)
    
    #priors
    n ~ dpois(lambda)
    p ~ dbeta(a,b)
    theta <- n*p
}")
inits <- list(p = 0.8,n = 1)
model <- jags.model(model_string,data = data,inits = inits, quiet = TRUE)
# burn in
update(model,10000, progress.bar = "none")
params <- c("n","p","theta")
samples <- coda.samples(model = model,variable.names = params,n.iter = 20000,progress.bar = "none")

summary(samples)
plot(samples)


### partc
Y <- 10
a <- 10
b <- 10
lambda <- 10

data <- list(Y = Y, a = a,b = b, lambda = lambda)
model_string <- textConnection("model{
    # liklihood
    Y ~ dbin(p,n)
    
    #priors
    n ~ dpois(lambda)
    p ~ dbeta(a,b)
    theta <- n*p
}")

inits <- list(p = 0.8,n = 1)
model <- jags.model(model_string,data = data,inits = inits, quiet = TRUE)
# burn in
update(model,10000, progress.bar = "none")
params <- c("n","p","theta")
samples <- coda.samples(model = model,variable.names = params,n.iter = 20000,progress.bar = "none")

summary(samples)
plot(samples)