# simulation setup from Sherwood and Wang
library(MASS)
library(penplaqr)
library(doParallel)
library(doRNG)

# setting sample size and dimension of predictors
# note that p > n in some cases


n <- 300



lambda_vals <- seq(from = 20, to = 70, length.out = 100)
set.seed(30)
beta_vals <- runif(4, 0.5, 1.5)

p <- 600

row_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2)
col_mat <- matrix(1:(p+2), nrow = p+2, ncol = p+2, byrow = T)

index_mat <- row_mat - col_mat

index_mat <- abs(index_mat)

big_sig <- 0.5^index_mat

diag(big_sig) <- 1

registerDoParallel(cores = 4)
registerDoRNG(43)
B <- 1:100

combos <- expand.grid(B, p_vec)

sim_results <- foreach(sim = B) %dopar%{


X <- mvrnorm(n, mu = rep(0, p+2), Sigma = big_sig)



# simulating covariate data from a multivariate normal distribution





errors_norm <- rnorm(n)





errors_tdist <- rt(n, df = 3)





psi <- rnorm(n, mean = 0, sd = 0.7)

errors_normhetero <- X[,1]*psi



# transforming covariates to make non-linear components
X[,1] <- sqrt(12)*pnorm(X[,1])
X[,25] <- pnorm(X[,25])
X[,26] <- pnorm(X[,26])


y_norm <- X[,c(6, 12, 15, 20)]%*%beta_vals + sin(2*pi*X[,25]) + X[,26]^3 + errors_norm

y_tdist <- X[,c(6, 12, 15, 20)]%*%beta_vals + sin(2*pi*X[,25]) + X[,26]^3 + errors_tdist

y_normhetero <- X[,c(6, 12, 15, 20)]%*%beta_vals + sin(2*pi*X[,25]) + X[,26]^3 + errors_normhetero

my_dat_norm <- as.data.frame(cbind(y_norm, X))

my_dat_tdist <- as.data.frame(cbind(y_tdist, X))

my_dat_normhetero <- as.data.frame(cbind(y_normhetero, X))

names(my_dat_norm) <- names(my_dat_tdist) <- names(my_dat_normhetero) <- c("y", paste0("V", 1:(p + 2)))

lin_formula <- as.formula(paste("y ~", paste0("V", c(1:24, 27:(p+2)), collapse = " + ")))
nonlin_formula <- as.formula(paste("~", paste0("V", c(25, 26), collapse = " + ")))



coef_array <- array(NA, dim = c(length(lambda_vals), p + 7, 2))
qbic_mat <- matrix(NA, nrow = length(lambda_vals), ncol = 2)

for(lambda in seq_along(lambda_vals)){



      model_norm <- try(penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_norm, tau = .5, lambda = lambda, a= 3.7,
                          penalty = "SCAD", method = "fn"))
      model_tdist <- try(penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_tdist, tau = .5, lambda = lambda, a = 3.7,
                              penalty = "SCAD", method = "fn"))

     # model_normhetero <- penplaqr(lin_formula, nonlinVars = nonlin_formula, data = my_dat_normhetero, tau = .7, lambda = lambda_vals[i],
                                  # a = 3.7, penalty = "SCAD", method = "br")
      if(class(model_norm) != "try-error"){
        coef_array[lambda, , 1] <- coef(model_norm)

        resids <- residuals(model_norm)[1:300]

        deg_free <- sum(abs(resids) < 1e-06)

        qbic_mat[lambda, 1] <- QBIC(.5, resids, p, deg_free)

      }


      if(class(model_tdist) != "try-error"){

        coef_array[lambda, , 2] <- coef(model_tdist)
        resids <- residuals(model_tdist)[1:300]
        deg_free <- sum(abs(resids) < 1e-06)
        print(deg_free)
        qbic_mat[lambda, 2] <- QBIC(.5, resids, p, deg_free)


      }

      }

  cbind(coef_array[qbic_mat[,1] == min(qbic_mat[,1]), , 1],
        coef_array[qbic_mat[,2] == min(qbic_mat[,2]), , 2])

}

true_indices <- c(7, 13, 16, 21)
beta0 <- rep(0, 301)
beta0[true_indices] <- beta_vals
mean(unlist(lapply(sim_results, function(x) sum(x[true_indices, 1] > 1e-8))))

mean(unlist(lapply(sim_results, function(x) sum(x[-c(1, true_indices, 302:307), 1] > 1e-8))))

mean(unlist(lapply(sim_results, function(x) sum(x[true_indices, 1] > 1e-8) == 4 & sum(x[-c(1, true_indices, 302:307), 1] > 1e-8) == 0)))

mean(unlist(lapply(sim_results, function(x) sum(x[true_indices, 2] > 1e-8))))

mean(unlist(lapply(sim_results, function(x) sum(x[-c(1, true_indices, 302:307), 2] > 1e-8))))

mean(unlist(lapply(sim_results, function(x) sum(x[true_indices, 2] > 1e-8) == 4 & sum(x[-c(1, true_indices, 302:307), 2] > 1e-8) == 0)))

mean(unlist(lapply(sim_results, function(x) sum((x[1:301, 1] - beta0)^2))))

mean(unlist(lapply(sim_results, function(x) sum((x[1:301, 2] - beta0)^2))))

mean(unlist(lapply(sim_results, function(x) x[2, 1] > 1e-08)))

mean(unlist(lapply(sim_results, function(x) x[2,2] > 1e-08)))

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




