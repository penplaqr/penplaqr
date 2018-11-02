# simulation setup from Sherwood and Wang
library(MASS)
library(penplaqr)



# setting sample size and dimension of predictors
# note that p > n
B <- 100
p_vec <- c(100, 300, 600)
n <- 300



lambda_vals <- seq(from = 10, to = 55, length.out = 100)

set.seed(30)
beta_vals <- runif(4, 0.5, 1.5)
error_dist <- c("normal_constantvar", "normal_nonconstantvar", "t_dist")

coef_array <- array(NA, dim = c(B, max(p_vec), length(p_vec), length(lambda_vals), length(error_dist)))


for(i in 1:B){
  for(p_val in seq_along(p_vec)){
    for(j in seq_along(error_dist)){
# simulating covariate data from a multivariate normal distribution

  p <- p_vec[p_val]

  row_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2)
  col_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2, byrow = T)

  index_mat <- row_mat - col_mat

  index_mat <- abs(index_mat)

  big_sig <- 0.5^index_mat

  diag(big_sig) <- 1

  X <- mvrnorm(n, mu = rep(0, p+2), Sigma = big_sig)


# simulating outcome data
# the first 150 covariates are unrelated to the outcome
# the next 45 covariates are linearly related to the outcome
# the final 5 covariates are related to the outcome through nonlinear functions

  if(error_dist[j] == "normal_constantvar"){

    errors <- rnorm(n)

  }

  if(error_dist[j] == "t_dist"){

    errors <- rt(n, df = 3)

  }

  if(error_dist[j] == "normal_nonconstantvar"){

    psi <- rnorm(n, mean = 0, sd = 0.7)

    errors <- X[,1]*psi

  }

# transforming covariates to make non-linear components
    X[,1] <- sqrt(12)*pnorm(X[,1])
    X[,25] <- sin(2*pi*pnorm(X[,25]))
    X[,26] <- pnorm(X[,26])^3


    y <- X[,c(6, 12, 15, 20)]%*%beta_vals + X[, 25:26]%*%c(1,1) + errors

    my_dat <- as.data.frame(cbind(y, X))

    names(my_dat) <- c("y", paste0("V", 1:(p + 2)))

    lin_formula <- as.formula(paste("y ~", paste0("V", c(1:24, 27:(p+2)), collapse = " + ")))
    nonlin_formula <- as.formula(paste("~", paste0("V", c(25, 26), collapse = " + ")))
    for(k in seq_along(lambda_vals)){
       model <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat, tau = ifelse(error_dist[j] == "normal_nonconstantvar", .7, .5), lambda = lambda_vals[k], a= 3.7,
                              penalty = "SCAD", method = "br")

        coef_array[i,1:p, p_val ,k,j] <- coef(model)[paste0("V", c(1:24, 27:(p+2)))]
        print(paste("did iteration", i, "at lambda value", k, "for error distribution", error_dist[j], "and p equal to", p))
    }
  }
  }
}

threshold <- 1e-5
zero_coefs <- lapply(model_list, function(x) sum(abs(coef(x)) < threshold))

mean(apply(coef_array[,,1], 1, function(x) sum((c(rep(0, 190), beta_vals) - x)^2)))


plot(x = lambda_vals, y = zero_coefs, type = 'b')
abline(h = 150)
abline(h = 195, col = 'red')
