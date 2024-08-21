#### assgn-5 que1
library(geoR)
library(rjags)
data("gambia")

gam <- gambia

Y <- gambia$pos
X <- as.matrix(gambia[,-3])
X <- scale(X)



data <- list(Y =Y,X = X, n = length(Y),p = dim(X)[2])
model_string <- textConnection("model{
# Likelihood 
for(i in 1:n){
   Y[i] ~ dbern(q[i])
   logit(q[i]) <- inprod(X[i,],beta[])
   }

#prior
for(j in 1:p){
     beta[j] ~ dnorm(0,0.01)
   }
}")
model <- jags.model(model_string, data = data, n.chains = 2,quiet = TRUE)
update(model,1000,progress.bar = "none")
param <- c("beta")
samples <- coda.samples(model, variable.names = param, 
                        n.iter = 5000,progress.bar = "none")

save(samples,file = "samples1.RData")

summary(samples)
load("samples1.RData")
plot(samples)




###que 1 part b
a <- 0.1
b <- 0.1
L <- 65
data("gambia")

gam <- gambia

Y <- gambia$pos
X <- as.matrix(gambia[,-3])
X <- scale(X)

ar <- as.data.frame(gambia$x)
cont <- table(ar)
s <- rep(1:65, times = c(cont))

data <- list(Y =Y,X = X,s = s, n = length(Y),p = dim(X)[2], L = L, a= a, b= b)
model_string <- textConnection("model{
# Likelihood 
for(i in 1:n){
   Y[i] ~ dbern(q[i])
   
   logit(q[i]) <- inprod(X[i,],beta[]) + alpha[s[i]]
   }

#prior
for(j in 1:p){
     beta[j] ~ dnorm(0,0.01)
  }
for(j in 1 : L){
    alpha[j] ~ dnorm(0,tau)
}   
 
 sigma <- 1/sqrt(tau)
 tau ~ dgamma(a,b) 
   
}")

model <- jags.model(model_string, data = data, n.chains = 2,quiet = TRUE)
update(model,1000,progress.bar = "none")
param <- c("beta")
samples <- coda.samples(model, variable.names = param, 
                        n.iter = 5000,progress.bar = "none")

save(samples,file = "samples2.RData")

summary(samples)

plot(samples)
load("samples2.RData")
plot(gambia$x, gambia$y)


### Que2

library(MASS)
library(rjags)

# Define the JAGS model
model_string <-textConnection( "
model {
  # Likelihood
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[z[i]], tau[z[i]])
    z[i] ~ dcat(p)
  }
    #priors
    for (k in 1:K) {
    mu[k] ~ dnorm(0, 0.0001)
    tau[k] ~ dgamma(0.01, 0.01)
    
    }
    
    # Mixing proportions with uniform Dirichlet prior
  
  
    
    p[1:K] ~ ddirch(alpha[1:K])
    
  
}")

# Prepare the data
Y <- galaxies
N <- length(Y)




K <- 3
alpha <- rep(1, K) 



# Run the model
data <- list(Y = Y, N = N, K = K, alpha = alpha)
model <- jags.model(model_string, data = data, 
                    n.chains = 2)


update(model, n.iter = 1000)
samples <- coda.samples(model, variable.names = c("mu","tau","p"),
                        n.iter = 5000, thin = 5)
summary(samples)



# Extract posterior samples
posterior_samples <- as.matrix(samples)
mu <- as.matrix(posterior_samples[,1:3])
tau <- as.matrix(posterior_samples[,7:9])

p <- as.matrix(posterior_samples[,4:6])

S <- 5000
val <- runif(n = S*351, min = 5000,max = 40000)
y <- matrix(val,nrow = S,ncol = 351)


den <- matrix(0,nrow = S,ncol = 351) 

for (i in 1:S) {
  for (j in 1:351) {
    den[i,j] <- p[i,1] *dnorm(y[i,j],mean = mu[i,1],sd = 1/tau[i,1]^2) + 
      p[i,2] *dnorm(y[i,j],mean = mu[i,2],sd = 1/tau[i,2]^2) + 
      p[i,3] *dnorm(y[i,j],mean = mu[i,3],sd = 1/tau[i,3]^2) 
  }
  
}

### que3
library(rjags)
Y <- c(563,10)
N_s <- c(2820,27)
c <- 1
n <- 2





# Baye factor
lgconst <- lgamma(Y[1] + Y[2] + 1) + (Y[1] + 1) * log(N_s[1]) + (Y[2] + 1) * log(N_s[2]) -
  ((Y[1] + Y[2] + 1) * log(N_s[1] + N_s[2]) + lgamma(Y[1] + 1) + lgamma(Y[2] + 1))

# Calculate the exponential of lgconst
result <- exp(lgconst)

bf21_c1 <- result*c*pgamma(c,shape = Y[1]+Y[2]+1, rate = N_s[1] + N_s[2])/(pgamma(c,Y[1] + 1,N_s[1])*pgamma(c,Y[2]+1,N_s[2]) )



data <- list(Y = Y, n = length(Y), N_s = N_s,c= c)
model_string <- textConnection("model{
# Likelihood 
for(i in 1:n){
Y[i] ~ dpois(N_s[i]*lambda[i])
like[i] <- dpois(Y[i],N_s[i]*lambda[i])
}

#prior
for(i in 1:n){
  lambda[i] ~ dunif(0,c)
}
}")


model <- jags.model(model_string, data = data, n.chains = 2,quiet = TRUE)
update(model,10000,progress.bar = "none")

samples <- coda.samples(model, variable.names = c("like"), 
                        n.iter = 50000,progress.bar = "none")

## dic
DIC1 <- dic.samples(model,n.iter = 50000 ,progress.bar = "none")
DIC1
# Compute WAIC

like <- rbind(samples[[1]],samples[[2]]) # Combine the two chains
fbar <- colMeans(like)
Pw <- sum(apply(log(like),2,var))
WAIC1 <- -2*sum(log(fbar)) + 2*Pw
WAIC1



### c=10
c <- 10
# Baye factor
lgconst <- lgamma(Y[1] + Y[2] + 1) + (Y[1] + 1) * log(N_s[1]) + (Y[2] + 1) * log(N_s[2]) -
  ((Y[1] + Y[2] + 1) * log(N_s[1] + N_s[2]) + lgamma(Y[1] + 1) + lgamma(Y[2] + 1))

# Calculate the exponential of lgconst
result <- exp(lgconst)

bf21_c10 <- result * c * pgamma(c,shape = Y[1]+Y[2]+1, rate = N_s[1] + N_s[2])/(pgamma(c,Y[1] + 1,N_s[1])*pgamma(c,Y[2]+1,N_s[2]) )





data <- list(Y = Y, n = length(Y), N_s = N_s,c= c)
model_string <- textConnection("model{
# Likelihood 
for(i in 1:n){
Y[i] ~ dpois(N_s[i]*lambda[i])
like[i] <- dpois(Y[i],N_s[i]*lambda[i])
}

#prior
for(i in 1:n){
  lambda[i] ~ dunif(0,c)
}
}")

model <- jags.model(model_string, data = data, n.chains = 2,quiet = TRUE)
update(model,10000,progress.bar = "none")

samples_10 <- coda.samples(model, variable.names = c("like"), 
                           n.iter = 50000,progress.bar = "none")

## dic
DIC_10 <- dic.samples(model,n.iter = 50000 ,progress.bar = "none")
DIC_10
# Compute WAIC

like <- rbind(samples_10[[1]],samples_10[[2]]) # Combine the two chains
fbar <- colMeans(like)
Pw <- sum(apply(log(like),2,var))
WAIC_10 <- -2*sum(log(fbar)) + 2*Pw
WAIC_10

BF21 <- c(bf21_c1,bf21_c10)
C <- c(1,10)
DIC <- c(14.35,14.35)
Waic <- c(WAIC1,WAIC_10)
tab <- data.frame(C,BF21,DIC,Waic)
print(tab)





### Assign-5 que-4


library(geoR)
library(rjags)
Y <- gambia$pos


X   <- gambia[,4:8]
X   <- scale(X)
data  <- list(Y=Y,X=X,n=length(Y))


# Fit logistic model
model_string <- textConnection("model{
 for(i in 1:n){
  Y[i] ~ dbern(pi[i])
  logit(pi[i]) <- beta[1] + X[i,1]*beta[2] +
                  X[i,2]*beta[3] + X[i,3]*beta[4] +
                  X[i,4]*beta[5] + X[i,5]*beta[6]
  
 }
 
 #priors
 for(j in 1:6){
 beta[j] ~ dnorm(0,0.01)
 }
 
 # Posterior preditive checks
for(i in 1:n){
Y1[i] ~ dbern(pi[i])
}
D[1] <- mean(Y1[]^2) - mean(Y1[])^2
D[2] <- mean(Y1[])

 
}")


model <- jags.model(model_string,data = data, n.chains=1,quiet=TRUE)

update(model, 500, progress.bar = "none")

samples <- coda.samples(model, variable.names = c("D"),
                       n.iter = 5000, thin = 5, progress.bar = "none")
plot(samples)

D.m <- samples[[1]]

##Compute the Bayesian p-values
D0 <- c(var(Y), mean(Y))
Dnames <- c("Var Y", "Mean Y")


# Compute the test stats for the models
pval1 <- rep(NA, 2)
names(pval1) <-  c("Var Y", "Mean Y")



for(j in 1:2){

 pval1[j] <- mean(D.m[ , j] > D0[j])
  
}

print(pval1)



### assign5 ques5
Y <- WWWusage
Y <- as.numeric(Y)
data1 <- list( Y = Y ) 


## ar1
Ar1_model <- textConnection("model{
 for(t in 2:100){
 Y[t] ~ dnorm(mu[t],tau)
 mu[t] <- beta[1] + beta[2]*Y[t-1]
 like[t] <- dnorm(Y[t],mu[t],tau)
 }
 ## priors
 ## priors
 for(j in 1:2){
 beta[j] ~ dnorm(0,0.01)
 }
 
 tau ~ dgamma(0.1,0.1)
 
 

 
 
}")


model <- jags.model(Ar1_model,data = data1, n.chains=2,quiet=TRUE)

update(model, 10000, progress.bar = "none")

## dic
dic1 <- dic.samples(model,n.iter = 50000)
print(dic1)

samps1 <- coda.samples(model, variable.names = c("like"),
                        n.iter = 50000, thin = 5, progress.bar = "none")

# Compute WAIC
 
 like <- rbind(samps1[[1]],samps1[[2]]) # Combine the two chains
 fbar <- colMeans(like)
 Pw <- sum(apply(log(like),2,var))
 WAIC1 <- -2*sum(log(fbar)) + 2*Pw
WAIC1




## ar2
Ar2_model <- textConnection("model{
 for(t in 3:100){
 Y[t] ~ dnorm(mu[t],tau)
 mu[t] <- beta[1] + beta[2]*Y[t-1] + beta[3]*Y[t-2]
 like[t] <- dnorm(Y[t],mu[t],tau)
 }
 ## priors
 for(j in 1:3){
 beta[j] ~ dnorm(0,0.01)
 }
 
 tau ~ dgamma(0.1,0.1)
 
 
 
 
}")

model <- jags.model(Ar2_model,data = data1, n.chains=2,quiet=TRUE)

update(model, 10000, progress.bar = "none")

## dic
dic2 <- dic.samples(model,n.iter = 50000)
print(dic2)


# Compute WAIC
samps2 <- coda.samples(model, variable.names = c("like"),
                       n.iter = 50000, thin = 5, progress.bar = "none")

like <- rbind(samps2[[1]],samps2[[2]]) # Combine the two chains
fbar <- colMeans(like)
Pw <- sum(apply(log(like),2,var))
WAIC2 <- -2*sum(log(fbar)) + 2*Pw
WAIC2



## ar3
Ar3_model <- textConnection("model{
 for(t in 4:100){
 Y[t] ~ dnorm(mu[t],tau)
 mu[t] <- beta[1] + beta[2]*Y[t-1] + beta[3]*Y[t-2] + beta[4]*Y[t-3]
 like[t] <- dnorm(Y[t],mu[t],tau)
 }
 ## priors
 for(j in 1:4){
 beta[j] ~ dnorm(0,0.01)
 }
 
 tau ~ dgamma(0.1,0.1)
 
 
}")

model <- jags.model(Ar3_model,data = data1, n.chains=2,quiet=TRUE)

update(model, 10000, progress.bar = "none")

## dic
dic3 <- dic.samples(model,n.iter = 50000)
print(dic3)


# Compute WAIC
samps3 <- coda.samples(model, variable.names = c("like"),
                       n.iter = 50000, thin = 5, progress.bar = "none")

like <- rbind(samps3[[1]],samps3[[2]]) # Combine the two chains
fbar <- colMeans(like)
Pw <- sum(apply(log(like),2,var))
WAIC3 <- -2*sum(log(fbar)) + 2*Pw
WAIC3



## ar4
Ar4_model <- textConnection("model{
 for(t in 5:100){
 Y[t] ~ dnorm(mu[t],tau)
 mu[t] <- beta[1] + beta[2]*Y[t-1] + beta[3]*Y[t-2] + beta[4]*Y[t-3] + beta[5]*Y[t-4]
 like[t] <- dnorm(Y[t],mu[t],tau)
 
 }
 ## priors
 for(j in 1:5){
 beta[j] ~ dnorm(0,0.01)
 }
 
 tau ~ dgamma(0.1,0.1)
 
}")

model <- jags.model(Ar4_model,data = data1, n.chains=2,quiet=TRUE)

update(model, 10000, progress.bar = "none")

## dic
dic4 <- dic.samples(model,n.iter = 50000)
print(dic4)


# Compute WAIC
samps4 <- coda.samples(model, variable.names = c("like"),
                       n.iter = 50000, thin = 5, progress.bar = "none")

like <- rbind(samps4[[1]],samps4[[2]]) # Combine the two chains
fbar <- colMeans(like)
Pw <- sum(apply(log(like),2,var))
WAIC4 <- -2*sum(log(fbar)) + 2*Pw
WAIC4


### table for comparision
L <- c(1:4)
DIC <- c(626.5,517.9 ,506.9,487.8)
WAIC <- c(WAIC1,WAIC2,WAIC3,WAIC4)
model_compar <- data.frame(L,DIC,WAIC)
print(model_compar)
