library(profvis)
library(lineprof)

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


my_dat <- as.data.frame(cbind(y, X))

names(my_dat) <- c("y", paste0("V", 1:200))

lin_formula <- as.formula(paste("y ~", paste0("V", 1:195, collapse = " + ")))
nonlin_formula <- as.formula(paste("~", paste0("V", 196:p, collapse = " + ")))



prof <- lineprof(pen_fit <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat, tau = .5, lambda = 12, a= 3.7,
                    penalty = "SCAD", method = "fn"))
