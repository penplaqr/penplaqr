

#' Check function.
#'
#' \code{check} computes the "check" loss function used in quantile regression.
#'
#' @param tau a numeric-valued scalar between 0 and 1.
#' @param x a numeric-valued scalar or vector.
#'
#' @return a numeric or vector of numerics of length \code{length(tau)}.
check <- function(tau, x){

  x*(tau - (x < 0))



}


#' High-dimensional BIC.
#'
#' \code{QBIC} computes the high-dimensional BIC described in Sherwood and Wang (2016).
#'
#' @param lambda a nonnegative scalar; the tuning parameter.
#' @param resids a vector of residuals from the candidate model.
#' @param p_linear the number of candidate linear covariates.
#' @param degree the degrees of freedom of the candidate model.
QBIC <- function(lambda, resids, p_linear, degfree){

  n <- length(resids)

  log(sum(check(resids))) - degfree*log(p_linear)*log(log(n))/(2*n)

}


#' Augmented Data.
#'
#' \code{augment_data} creates an augmented dataset for use in the weighted
#' unpenalized quantile regression as part of the algorithm described in
#' Sherwood and Wang (2016) for solving the penalized quantile regression problem.
#'
#' @param n_linvars the number of predictors modeled with linear terms
#' @param n_nonlinvars the number of predictors modeled non-linearly
augment_data <- function(n_linvars, n_nonlinvars){

  extraobs_firstset <- extraobs_secondset <-  matrix(0, nrow = n_linvars, ncol = (1 + n_linvars + n_nonlinvars))
  extraobs_firstset[, 2:(n_linvars + 1)] <- diag(1, n_linvars, n_linvars)

  extraobs_secondset[, 2:(n_linvars + 1)] <- diag(-1, n_linvars, n_linvars)

  as.data.frame(rbind(extraobs_firstset, extraobs_secondset))



}


#' Penalty weights.
#'
#' \code{penalty_weights} calculates the weights associated with the
#' observations in the augmented dataset used in fitting the penalized partially
#' linear quantile regression.
#' @param beta_vec a numeric vector containing the current estimates of the
#' regression parameters from the linear portion of the model.
#' @param n the number of observations in the data set.
#' @param penalty_type a character string indicating which penalty function
#' should be used. See Details for more information.
#' @param penalty_deriv a user-specified function for use in creating the weights.
#' See Details for more information.
#' @param lambda a scalar corresponding to the tuning parameter.
#' @param a a scalar for use in the SCAD and MCP penalty functions.
#' @param ... (optional) additional arguments to be passed to the user-supplied penalty derivative function
#'
#' @details The weights are based on the derivative of the penalty function
#' used in defining the penalized regression problem. There are two ways for specifying
#' the derivative of the penalty function. First, the user can provide the name of the penalty
#' function through the argument penalty_type. Based off this specification, the derivative is then
#' computed using built-in functions. The current options for the penalty function are
#' "SCAD" and "MCP". Alternatively, the user can provide their own
#' penalty derivative function. This function should take as an argument a numeric vector
#' and return a numeric vector of the same length.
#'
#'
#' @return a vector of length (n + 2*length(beta_vec)) containing the weights associated
#' with the observations in the augmented data set.
penalty_weights <- function(beta_vec, n,  penalty_type = NULL, penalty_deriv = NULL, lambda, a = NULL, ... ){

  if(is.null(penalty_type) & is.null(penalty_deriv)) stop("You must either specify a penalty function or provide a penalty derivative")

  if(!(penalty_type %in% c("LASSO", "SCAD", "MCP"))) stop("The penalty function you specified is not supported")

  if(penalty_type == "MCP")  deriv_vals <- rqPen::mcp_deriv(beta_vec, lambda, a)

  if(penalty_type == "SCAD") deriv_vals <- rqPen::scad_deriv(beta_vec, lambda, a)

  if(!is.null(penalty_deriv)) deriv_vals <- penalty_deriv(beta_vec, lambda, ...)


  c(rep(1, n), deriv_vals, deriv_vals)

}

#' Partially linear penalized quantile regression
#'
#' \code{penplaqr} fits a partially linear penalized quantile regression
#' @param formula a formula object containing the response variable on the LHS and the variables to be modeled with
#' linear terms on the RHS.
#' @param nonlinVars a formula object with an empty LHS and the variables to be modeled non-linearly on the RHS.
#' @param tau a scalar indicating the desired quantile.
#' @param lambda a scalar indicating the value of the tuning parameter.
#' @param penalty a character string indicating which penalty should be used in fitting the model.
#' @param penalty_deriv an optional user-specified function for computing the derivative of the penalty function.
#' @param a additional tuning paramater for certain penalty functions.
#' @param epsilon scalar value for determining when convergence has been reached.
#' @param data a data frame
#' @import quantreg
#' @import splines
#' @export
penplaqr <- function(formula, nonlinVars = NULL, tau = .5, lambda = NULL, penalty = "SCAD", penalty_deriv = NULL, a = NULL, data = NULL, subset,
                      na.action, method = c("br", "fn"), model = TRUE, contrasts = NULL, init_beta = NULL,
                     splinesettings = NULL, epsilon = 1e-6, ...){
  if(length(tau)>1) stop("tau must be a number (not a vector) strictly between 0 and 1.")
  if(class(nonlinVars)!="formula" & !is.null(nonlinVars)){
    stop("nonlinVars must be of class \"formula\" or NULL. NULL by default.\n")
  }
  penplaqrcall <- match.call()
  rqcall <- penplaqrcall
  nonlinvars <- NULL
  linvars <- attr(terms(formula, data=data), "term.labels")

  # creating spline terms for non-linear portion of the model
  if(is.null(nonlinVars)){
    int <- ifelse(attr(terms(formula, data=data),"intercept")==1, "1", "0")
    rqcall$formula <- update(formula, paste(c("~",linvars,int),
                                            collapse="+"))
  } else {
    nonlinvars <- attr(terms(nonlinVars, data=data), "term.labels")
    nonlinvars <- nonlinvars[!(nonlinvars %in% all.vars(formula)[1])]
    linvars <- linvars[!(linvars %in% nonlinvars)]
    nonlinvarsbs <- rep(NA, length(nonlinvars))
    if(length(splinesettings)!=length(nonlinvarsbs) && !is.null(splinesettings))
      stop("splinesettings must either be NULL or a list with length = number of nonlinear covariates \n (see details ?plaqr)")

    bslist <- vector("list", length(nonlinvars))
    for( i in 1:length(nonlinvars) ){
      bslist[[i]] <- paste("bs(",nonlinvars[i], sep="")

      if( !is.null(splinesettings[[i]]$df) )
        bslist[[i]] <- paste(bslist[[i]], ", df=",splinesettings[[i]]$df, sep="")
      if( !is.null(splinesettings[[i]]$knots) )
        bslist[[i]] <- paste(bslist[[i]], ", knots=c(",
                             toString(splinesettings[[i]]$knots), ")", sep="")
      if( !is.null(splinesettings[[i]]$degree) )
        bslist[[i]] <- paste(bslist[[i]], ", degree=",splinesettings[[i]]$degree, sep="")
      if( !is.null(splinesettings[[i]]$Boundary.knots) )
        bslist[[i]] <- paste(bslist[[i]], ", Boundary.knots=c(",
                             toString(splinesettings[[i]]$Boundary.knots), ")", sep="")

      bslist[[i]] <- paste(bslist[[i]], ")", sep="")
    }

    nonlinvarsbs <- do.call(c, bslist)
    rqcall$formula <- update(formula, paste(c("~","1",linvars,nonlinvarsbs),
                                            collapse="+"))
  }

  # creating augmented dataset
  linear_terms <- attr(terms(formula, data = data), "term.labels")

  nonlinear_terms <- attr(terms(nonlinVars, data = data), "term.labels")

  response_var <- all.vars(update(formula, .~1))


    # removing variables included in non-linear terms from set of
    # linear terms
  linear_terms <- setdiff(linear_terms, nonlinear_terms)


  temp_dat <- data[,c(response_var, linear_terms, nonlinear_terms)]

  extra_observations <- augment_data(length(linear_terms), length(nonlinear_terms))

  names(extra_observations) <- names(temp_dat)


  augmented_data <- rbind(temp_dat, extra_observations)

  # calculating weights at initial value for beta
  if(is.null(init_beta)) init_beta <- rep(0, length(linear_terms))

  init_weights <- penalty_weights(abs(init_beta), n = nrow(data),
                                  penalty_type = penalty, penalty_deriv = penalty_deriv, lambda = lambda, a = a)

  # adding augmented data and weights to rq call
  rqcall[[1]] <- as.name("rq")
  rqcall$nonlinVars <- rqcall$splinesettings <- rqcall$lambda <- rqcall$penalty <- rqcall$epsilon <-  NULL
  rqcall$data <- augmented_data
  rqcall$method <- method
  rqcall$weights <- init_weights
  rqcall$a <- NULL
  # running first iteration
  first_fit <- eval.parent(rqcall)
  class(first_fit) <- "rq"

  current_beta <- summary(first_fit)$coefficients[linear_terms, 1]
  current_beta <- coef(first_fit)[linear_terms]
  prev_beta <- init_beta
  iter <- 1

  # iterating until convergence
  while(sum((current_beta - prev_beta)^2) > epsilon){

   prev_beta <- current_beta

   current_weights <- penalty_weights(abs(prev_beta), n = nrow(data),
                                      penalty_type = penalty, penalty_deriv = penalty_deriv, lambda = lambda, a = a)

   rqcall$weights <- current_weights
   current_fit <- eval.parent(rqcall)
   class(current_fit) <- "rq"


   current_beta <- coef(current_fit)[linear_terms]

    iter <- iter + 1

  }

  if(iter == 1) return(first_fit)
  if(iter > 1) return(current_fit)

}




