#
# Question 1
#
# Use Newton’s method to find x ∈ (0, π/2) such that
# x − sin(x)/cos(x) = −x
# What is the significance of this particular value of x?
#

f <- function(x)
  2*x - sin(x)/cos(x)

df <- function(x)
  2 - 2/(cos(2*x) + 1)

newton <- function(x, f, df)
  x <- x - f(x)/df(x)

# For 0 or π/2 it does not converge
x <- pi/2
for(i in 1:25) {
  x <- newton(x, f, df)
  print(x)
}

#
# Question 2
#
# Solve the maximum expected returns optimization exercise in 
# section 7 of slides 8 (that is, reproduce the solution presented 
# in the lecture). The last line of your R program should display 
# the portfolio weights vector w.
#

G <- function(x, mu, Sigma, sigmaP2)
{
  # n = number of assets
  # x[1:5] = our n different weights
  n <- length(mu)
  c(mu + rep(x[n+1], n) + 2*x[n+2]*(Sigma %*% x[1:n]),
    sum(x[1:n]) - 1,
    t(x[1:n]) %*% Sigma %*% x[1:n] - sigmaP2)
}

DG <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  grad <- matrix(0.0, n+2, n+2)
  grad[1:n, 1:n] <- 2*x[n+2]*Sigma
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- 2*(Sigma %*% x[1:5])
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- 2*t(x[1:5]) %*% Sigma
  grad
}

x <- c(rep(0.5, 5), 1, 1)
mu <- c(0.08, 0.10, 0.13, 0.15, 0.20)
Sigma <- matrix(c(0.019600, -0.007560, 0.012880, 0.008750, -0.009800,
                  -0.007560, 0.032400, -0.004140, -0.009000, 0.009450,
                  0.012880, -0.004140, 0.052900, 0.020125, 0.020125,
                  0.008750, -0.009000, 0.020125, 0.062500, -0.013125,
                  -0.009800, 0.009450, 0.020125, -0.013125, 0.122500), 5, 5)

for(i in 1:25)
  x <- x - solve(DG(x, mu, Sigma, 0.25^2),
                 G(x, mu, Sigma, 0.25^2))

mu2 <- t(mu) %*% x[1:5]
# Weights of each stock
cat("Optimum weights: ", x[1:5])
cat("Expected return: ", mu2)

#
# Question 3
#
# Using the same expected returns vector μ and covariance matrix 
# Σ from the previous exercise, solve the optimization problem
#
# minimize: w^TΣw 
# subject to: e^Tw = 1
# μ^Tw = μ2
# 
# where μ2 is the expected return computed in the previous exercise. 
# Again, the last line of your R program should display the portfolio 
# weights vector w. Compare with the previous answer.
#

min_return <- function(mu, Sigma, mu2) {
  n = length(mu)
  df <- matrix(0.0, n+2, n+2)
  df[1:n, 1:n] <- 2*Sigma
  df[n+1, 1:n] <- rep(1, n)
  df[n+2, 1:n] <- mu
  df[1:n, 6] <- rep(1, n)
  df[1:n, 7] <- mu
  
  b <- c(rep(0, n), 1, mu2)
  
  solve(df, b)
}

solution3 <- min_return(mu, Sigma, mu2)
std_dev <- sqrt(t(solution3[1:5]) %*% Sigma %*% solution3[1:5])

cat("Optimum weights: ", solution3[1:5])
cat("Standard deviation: ", std_dev)