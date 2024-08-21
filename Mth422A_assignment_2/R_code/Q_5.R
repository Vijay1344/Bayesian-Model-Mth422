library(MASS)
library(ggplot2)

# parameters
mu <- c(0, 0)          # means
sigma <- matrix(c(1, 1, 1, 4), nrow = 2)  # covariance matrix
B <- 10^5              # number of samples

# samples
samples <- mvrnorm(n = B, mu = mu, Sigma = sigma)

# Create data frame
df <- data.frame(X = samples[,1], Y = samples[,2])

#  scatter plot with 2D-kernel density heatmap
ggplot(df, aes(x = X, y = Y)) +
  geom_point(alpha = 0.2) +
  geom_density_2d() +
  labs(x = "X", y = "Y", title = "Bivariate Normal Distribution with Correlation 0.5") +
  theme_minimal()

# initial conditions and parameters
B0 <- 10^3
B <- 10^5   
rho <- 0.5  # Correlation
X <- numeric(B0 + B)
Y <- numeric(B0 + B)
X[1] <- 0  # Initial value of X
Y[1] <- 0  # Initial value of Y


# functions for conditional distributions
conditional.X.given.Y <- function(Y, rho, mu1, mu2, sigma1, sigma2) {
  mu.x.given.y <- mu1 + rho * (sigma1/sigma2) * (Y - mu2)
  sd.x.given.y <- sqrt(sigma1^2 * (1 - rho^2))
  
  return(rnorm(1, mean = mu.x.given.y, sd = sd.x.given.y))
}

conditional.Y.given.X <- function(X, rho, mu1, mu2, sigma1, sigma2) {
  mu.y.given.x <- mu2 + rho * (sigma2/sigma1) * (X - mu1)
  sd.y.given.x <- sqrt(sigma2^2 * (1 - rho^2))
  
  return(rnorm(1, mean = mu.y.given.x, sd = sd.y.given.x))
}

#  sampling
for (b in 2:(B0 + B))
{
  X[b] <- conditional.X.given.Y(Y[b-1], rho, 0, 0, 1, 2) #  mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 2
  Y[b] <- conditional.Y.given.X(X[b], rho, 0, 0, 1, 2) #  mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 2
}

# scatter plot with 2D-kernel density heatmap by removing B0 sample
df_markov <- data.frame(X = X[(B0 + 1):(B0 + B)], Y = Y[(B0 + 1):(B0 + B)])

ggplot(df_markov, aes(x = X, y = Y)) +
  geom_point(alpha = 0.2) +
  geom_density_2d() +
  labs(x = "X", y = "Y", title = "Markov Chain Sampling from Bivariate Normal Distribution") +
  theme_minimal()


##1-D kernel estimate
## for X componet
kde <- density(df_markov$X)
plot(kde)
# for Y component
kde.y <- density(df_markov$Y)
plot(kde.y)