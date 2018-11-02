library(profvis)
library(lineprof)

library(rqPen)
library(penplaqr)
library(MASS)

x <- rnorm(100)

y <- 1 + 5*x + rnorm(100)

my_lm <- lm(y ~ x)


my_lm


A <- matrix(2, nrow = 2, ncol = 2)

B <- matrix(c(3,2), nrow = 2, ncol = 2)

p <- 5
row_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2)
col_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2, byrow = T)

index_mat <- row_mat - col_mat

index_mat <- abs(index_mat)

big_sig <- 0.5^index_mat

diag(big_sig) <- 1

p <- 200
n <- 500

# simulating covariate data from a multivariate normal distribution
set.seed(30)
X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(1, nrow = p))

# simulating outcome data
# the first 150 covariates are unrelated to the outcome
# the next 45 covariates are linearly related to the outcome
# the final 5 covariates are related to the outcome through nonlinear functions

y <- X[,151:195]%*%rnorm(45, 3, sd = 1) + sin(X[,196:200])%*%rep(1, 5) + log(X[,196:200]^2)%*%rep(1, 5) + rnorm(n)

lambda_vals <- seq(from = 1, to = 20, by = 1)
my_dat <- as.data.frame(cbind(y, X))

names(my_dat) <- c("y", paste0("V", 1:200))

lin_formula <- as.formula(paste("y ~", paste0("V", 1:195, collapse = " + ")))
nonlin_formula <- as.formula(paste("~", paste0("V", 196:p, collapse = " + ")))


tic1 <- Sys.time()
pen_fit <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat, tau = .5, lambda = 12, a= 3.7,
                  penalty = "SCAD", method = "fn")
toc1 <- Sys.time()


tic2 <- Sys.time()
pen_fit <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat, tau = .5, lambda = 12*n, a= 3.7,
                    penalty = "SCAD", method = "br")

toc2 <- Sys.time()

tic3 <- Sys.time()
rqpen_fit <- rq.nc.fit(x = X, y = y, lambda = 12)
toc3 <- Sys.time()


tic4 <- Sys.time()
rqpen_fit <- rq.nc.fit(x = X, y = y, lambda = 12, alg = "LP")
toc4 <- Sys.time()


getAnywhere(rq.nc.fit)
