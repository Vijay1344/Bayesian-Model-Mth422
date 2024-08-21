
### assgn4 - ques4
Y1 <- c(2,-3.1,-1,0.2,0.3,0.4)
Y2 <- c(-3.5,-1.6,-4.6,-0.9,-5.1,0.1)

# Statistics from group 1
Ybar1 <- mean(Y1)
s21 <- sum((Y1 - Ybar1)^2)/length(Y1)
n1 <- length(Y1)

# Statistics from group 2
Ybar2 <- mean(Y2)
s22 <- sum((Y2- Ybar2)^2)/length(Y2)
n2 <- length(Y2)

# Posterior of the difference assuming equal variance
delta_hat <- Ybar2 - Ybar1
delta_hat
s2 <- (n1*s21 + n2*s22)/(n1+n2)
scale <- sqrt(s2)*sqrt(1/n1 + 1/n2)
df <- n1+n2
cred_int <- delta_hat + scale * qt(c(0.025, 0.975), df = df)
cred_int

# Posterior of delta assuming unequal variance using MC sampling
mu1 <- Ybar1 + sqrt(s21/n1)*rt(1e5,df = df)
mu2 <- Ybar2 + sqrt(s22/n2)*rt(1e5,df = df)
delta <- mu2 - mu1

plot(density(delta), main = "Posterior distribution of the difference in means")
quantile(delta, c(0.025, 0.975)) # 95% credible set


#### sensitivity

data <- list(n = 6, Y1 = Y1, Y2 = Y2)
library(rjags)

model_string <- textConnection("model{
    # liklihood
    for(i in 1:n){
    Y1[i] ~ dnorm(mu, tau)
    Y2[i] ~ dnorm(mu + delta,tau)
    }
    
    #priors
     mu ~ dnorm(0,0.0001)
     delta ~ dnorm(0,0.0001)
     tau ~ dgamma(0.01,0.01)
}")


model <- jags.model(model_string,data = data, n.chains=2,quiet=TRUE)
update(model, 10000, progress.bar="none")
params  <- c("delta")
samples <- coda.samples(model, 
                        variable.names=params, 
                        n.iter=20000, progress.bar="none")



summary(samples)
plot(samples)