#' @title Count divergences
#'
#' @description This is a helper function to count the number of divergences
#' for an object of class [rstan::stanfit()]
#'
#' @param stan_fit an object of class [rstan::stanfit()]
#'
#' @return an integer giving the total number of divergences across all chains
#'
#' @references
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' @examples
#'
#' @importFrom rstan get_sampler_params

count_stan_divergences = function(stan_fit) {
  foo = get_sampler_params(stan_fit, inc_warmup = FALSE);
  n_draws = lapply(foo, nrow)[[1]];
  sum(unlist(lapply(foo,"[",i = 1:n_draws, j = "divergent__")));
}
