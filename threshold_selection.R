# simulation setup from Sherwood and Wang
library(MASS)
library(penplaqr)
library(doParallel)
library(doRNG)

# setting sample size and dimension of predictors
# note that p > n

p_vec <- c(100, 300, 600)
n <- 300



lambda_vals <- seq(from = 10, to = 55, by = 1)

set.seed(30)
beta_vals <- runif(4, 0.5, 1.5)

p <- 600
coef_array <- array(NA, dim = c(p + 2, length(lambda_vals), 2))

row_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2)
col_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2, byrow = T)

index_mat <- row_mat - col_mat

index_mat <- abs(index_mat)

big_sig <- 0.5^index_mat

diag(big_sig) <- 1

X <- mvrnorm(n, mu = rep(0, p+2), Sigma = big_sig)



# simulating covariate data from a multivariate normal distribution





errors_norm <- rnorm(n)





errors_tdist <- rt(n, df = 3)





psi <- rnorm(n, mean = 0, sd = 0.7)

errors_normhetero <- X[,1]*psi



# transforming covariates to make non-linear components
X[,1] <- sqrt(12)*pnorm(X[,1])
X[,25] <- sin(2*pi*pnorm(X[,25]))
X[,26] <- pnorm(X[,26])^3


y_norm <- X[,c(6, 12, 15, 20)]%*%beta_vals + X[, 25:26]%*%c(1,1) + errors_norm

y_tdist <- X[,c(6, 12, 15, 20)]%*%beta_vals + X[, 25:26]%*%c(1,1) + errors_tdist

y_normhetero <- X[,c(6, 12, 15, 20)]%*%beta_vals + X[, 25:26]%*%c(1,1) + errors_normhetero

my_dat_norm <- as.data.frame(cbind(y_norm, X))

my_dat_tdist <- as.data.frame(cbind(y_tdist, X))

my_dat_normhetero <- as.data.frame(cbind(y_normhetero, X))

names(my_dat_norm) <- names(my_dat_tdist) <- names(my_dat_normhetero) <- c("y", paste0("V", 1:(p + 2)))

lin_formula <- as.formula(paste("y ~", paste0("V", c(1:24, 27:(p+2)), collapse = " + ")))
nonlin_formula <- as.formula(paste("~", paste0("V", c(25, 26), collapse = " + ")))




registerDoParallel(cores = 4)
registerDoRNG(1235)


sims   <-   foreach(i = lambda_vals, .combine = rbind) %dopar% {

      # simulating outcome data
      # the first 150 covariates are unrelated to the outcome
      # the next 45 covariates are linearly related to the outcome
      # the final 5 covariates are related to the outcome through nonlinear functions



      model_norm <- try(penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_norm, tau = .5, lambda = i, a= 3.7,
                          penalty = "SCAD", method = "br"))
      model_tdist <- try(penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_tdist, tau = .5, lambda = i, a = 3.7,
                              penalty = "SCAD", method = "br"))

     # model_normhetero <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_normhetero, tau = .7, lambda = lambda_vals[i],
                                  # a = 3.7, penalty = "SCAD", method = "br")

      if(class(model_norm) == "try-error"){
        df_normal <- data.frame(as.list(c("normal", i, rep(NA, p))))
        names(df_normal) <- c("error", "lambda", paste0("V", (1:602)[-c(25, 26)]))
        df_normal[, grep("V", names(df_normal))] <- as.numeric(df_normal[, grep("V", names(df_normal))])
      }
      if(class(model_tdist) == "try-error"){
        df_tdist <- data.frame(as.list(c("tdist", i, rep(NA, p))))
        names(df_tdist) <- c("error", "lambda", paste0("V", (1:602)[-c(25, 26)]))
        df_tdist[, grep("V", names(df_tdist))] <- as.numeric(df_tdist[, grep("V", names(df_tdist))])
      }
      # coef_array[,i,3] <- coef(model_normhetero)[paste0("V", 1:102)]
      if(class(model_norm) != "try-error") df_normal <- cbind(data.frame(error = "normal", lambda = i ), as.data.frame( as.list(coef(model_norm)[paste0("V", (1:602)[-(25:26)])])))
      if(class(model_tdist) != "try-error") df_tdist <- cbind(data.frame(error = "tdist", lambda = i), as.data.frame(as.list(coef(model_tdist)[paste0("V", (1:602)[-(25:26)])])))

      rbind(df_normal, df_tdist)

      }

sims_normal <- sims[sims$error == "normal", ]
sims_tdist <- sims[sims$error == "tdist", ]
threshold <- 1e-08
p_nonzero_normal <- p - rowSums(sims_normal[, grep("V", names(sims_normal))] < threshold)
p_nonzero_tdist <- p - rowSums(sims_tdist[, grep("V", names(sims_normal))] < threshold)
# determining sparsity of each fit



sparsity_array <- coef_array
sparsity_array[abs(sparsity_array) < threshold] <- 0


cbind(apply(coef_array[,,1], 2, function(x) sum(abs(x) < threshold, na.rm = T)), lambda_vals)

rowSums(sims[sims$error == "tdist" & sims$lambda == 40, 3:ncol(sims)] < threshold)
rowSums(sims[sims$error == "tdist" & sims$lambda == 41, 3:ncol(sims)] < threshold)
n_zero <- apply(coef_array[,,1], 2, function(x) sum(abs(x) < threshold, na.rm = T))

summary(n_zero)




