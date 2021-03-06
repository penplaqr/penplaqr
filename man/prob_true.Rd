% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workHorse.R
\name{prob_true}
\alias{prob_true}
\title{Probability of Selecting the True Model.}
\usage{
prob_true(model_coefficients, true_coefficients, threshold)
}
\arguments{
\item{model_coefficients}{a vector or matrix of coefficient estimates. If a matrix,
then the coefficients should vary by row and the simulations should vary by column.}

\item{true_coefficients}{a vector of indices of the true non-zero coefficients.}

\item{threshold}{a scalar value below which a coefficient should be considered
equal to zero.}
}
\value{
a scalar indicating the proportion of simulations in which the true
model was selected exactly.
}
\description{
Given coefficient estimates from a set of simulations, \code{prob_true}
calculates the proportion of simulations that selected the true model exactly.
}
