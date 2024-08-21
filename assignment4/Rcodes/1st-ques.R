set.seed(2)

library(rjags)
Y <- c(64,13,33,18,30,20)
n <- length(Y)
data <- list(Y = Y, n=n)
model_string <- textConnection("model{
    # liklihood
    for(i in 1:n){
    Y[i] ~ dpois(exp(alpha + i*beta))
    }
    
    #priors
    alpha ~ dnorm(0,0.0001)
    beta ~ dnorm(0,0.0001)
}")

inits <- list( alpha = 0.5,beta = 2)
model <- jags.model(model_string,data = data,inits = inits,n.chains = 2,quiet = TRUE)

# burn in
update(model,10000, progress.bar = "none")
params <- c("alpha","beta")
samples <- coda.samples(model = model,variable.names = params,n.iter = 20000,progress.bar = "none")

summary(samples)
plot(samples)
