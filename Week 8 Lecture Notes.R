# Black-Scholes call price
bsc <- function(S, T, t, K, r, s, q) {
  d1 <- (log(S/K)+(r-q+0.5*s^2)*(T-t))/(s*sqrt(T-t))
  d2 <- d1-s*sqrt(T-t)
  S*exp(-q*(T-t))*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
}
# Lecture 8.1
bsc(50, 0.5, 0.0, 45, 0.06, 0.2, 0.02)

# Placing multiple values in for sigma
bsc(50, 0.5, 0.0, 45, 0.06, c(0.15, 0.2, 0.25), 0.02)

sigmas <- seq(0.05, 0.5, by = 0.01)
fsig <- bsc(50, 0.5, 0.0, 45, 0.06, sigmas, 0.02) - 7
plot(sigmas, fsig, type = "l")
bsc(50, 0.5, 0.0, 45, 0.06, 0.25, 0.02) - 7


# Lecture 8.2
bisection <- function(f, a, b, tol = 0.001) {
  while(b-a > tol) {
    c <- (a+b)/2
    if(sign(f(c)) == sign(f(a)))
      a <- c else
        b <- c }
  (a+b)/2 
}

fsig <- function(sigma)
  bsc(50, 0.5, 0.0, 45, 0.06, sigma, 0.02) - 7

bisection(fsig, 0.1, 0.3)

bsc(50, 0.5, 0.0, 45, 0.06, 0.2511719, 0.02)


# Lecture 8.4
F <- function(x, y)
  c(4*(x - 1)^3, 4*(y - 1)^3)

DF <- function(x, y)
  diag(c(12*(x - 1)^2, 12*(y - 1)^2))

x <- c(0, 0)

for(i in 1:25)
  x <- x - solve(DF(x[1], x[2]), F(x[1], x[2]))

# Lecture 8.5

F <- function(x)
  c(-4+2*x[4]*x[1],
    6+2*x[4]*x[2],
    2*x[1]-x[2]-x[3],
    x[1]^2+x[2]^2-13)

DF <- function(x)
  matrix(c(2*x[4],     0,   0, 2*x[1],
           0, 2*x[4],  0, 2*x[2],
           2,     -1, -1,      0,
           2*x[1], 2*x[2],  0,      0),
         4, 4, byrow = TRUE)

x <- rep(1, 4)
for(i in 1:15)
  x <- x - solve(DF(x), F(x))

x <- -rep(1, 4)
for(i in 1:15)
  x <- x - solve(DF(x), F(x))


# Lecture 8.6
G <- function(x)
  c(3 + 6*x[6]*x[1], -4 - 2*x[5]*x[2],
    1 + 2*x[5]*x[3] + 2*x[6]*x[3],
    -2 + 2*x[5]*x[4] + 4*x[6]*x[4],
    -x[2]^2 + x[3]^2 + x[4]^2 - 1,
    3*x[1]^2 + x[3]^2 + 2*x[4]^2 - 6)

DG <- function(x) {
  grad <- matrix(0, 6, 6)
  grad[1,] <- c(6*x[6], 0, 0, 0, 0, 6*x[1])
  grad[2,] <- c(0, -2*x[5], 0, 0, -2*x[2], 0)
  grad[3,] <- c(0, 0, 2*x[5] + 2*x[6], 0, 2*x[3], 2*x[3])
  grad[4,] <- c(0, 0, 0, 2*x[5] + 4*x[6], 2*x[4], 2*x[4])
  grad[5,] <- c(0, -2*x[2], 2*x[3], 2*x[4], 0, 0)
  grad[6,] <- c(6*x[1], 0, 2*x[3], 4*x[4], 0, 0)
  grad
}

x <- c(1, -1, 1, -1, 1, -1)

for(i in 1:25)
  x <- x - solve(DG(x), G(x))

f <- function(x) 3*x[1] - 4*x[2] + x[3] - 2*x[4]

f(x)

# Lecture 8.7

G <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  c(mu + rep(x[n+1], n) + 2*x[n+2]*(Sigma %*% x[1:n]),
    sum(x[1:n]) - 1,
    t(x[1:5]) %*% Sigma %*% x[1:5] - sigmaP2)
}

DG <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  grad <- matrix(0.0, n+2, n + 2)
  grad[1:n, 1:n] <- 2*x[n+2]*Sigma
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- 2*(Sigma %*% x[1:5])
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- 2*t(x[1:5]) %*% Sigma
  grad
}

x <- c(rep(0.5, 5), 1, 1)

for(i in 1:25)
  x <- x - solve(DG(x, mu, Sigma, 0.25^2),
                 G(x, mu, Sigma, 0.25^2))