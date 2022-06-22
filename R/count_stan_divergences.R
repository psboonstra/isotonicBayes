#' @title Count divergences
#'
#' @description This is a helper function to count the number of divergences
#' for an object of class [rstan::stanfit()]. It is used within `bayesian_isotonic`
#' but can also be used on its own
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
#' set.seed(1)
#' fake_data <- make_grouped_data(runif(100), rbinom(100, 1, 0.5))
#' fake_model <- bayesian_isotonic(data_grouped = fake_data, return_as_stan_object = TRUE)
#' count_stan_divergences(fake_model)
#'
#' @export
#'
#' @importFrom rstan get_sampler_params

count_stan_divergences = function(stan_fit) {
  foo = get_sampler_params(stan_fit, inc_warmup = FALSE);
  n_draws = lapply(foo, nrow)[[1]];
  sum(unlist(lapply(foo,"[",i = 1:n_draws, j = "divergent__")));
}
