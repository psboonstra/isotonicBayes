#' The 'isotonicBayes' package.
#'
#' @description This R package implements the isotonic probability vector priors
#' developed in Boonstra, Owen, and Kang (2020). The main function is
#' `bayesian_isotonic`, through which users can fit isotonic regression models
#' with a binary outcome regressed against a categorical predictor.
#'
#' @docType package
#' @name isotonicBayes-package
#' @aliases isotonicBayes
#' @useDynLib isotonicBayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#'
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
