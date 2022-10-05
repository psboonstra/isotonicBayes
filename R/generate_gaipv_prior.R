#' @title Generate draws from GAIPV prior
#'
#' @description Generate draws from the GAIPV prior and then calculate the
#' effective number of parameters according to this prior.
#'
#' @param n_cats_seq A vector of one or more positive integers giving the
#' desired values of K over which to evaluate
#'
#' @param shape_seq A vector of one or more positive numbers giving the
#' set of shape parameters s over which to evaluate
#'
#'
#' @param n_draws a positive integer giving the number of random draws from the
#' t-distributions with which to estimate expectations
#'
#' @param seed a positive integer giving the random seed, provided for reproducibility.
#'
#' @return a tibble with number of rows equal to the product of the lengths of
#' `n_cats_seq` and `shape_seq` and three columns: the current number of categories,
#' the current shape, and the calculated effective number of prior parameters
#'
#' @references
#'
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' \insertRef{piironen2017hyperprior}{isotonicBayes}
#'
#' @examples
#' generate_gaipv_prior(n_cats_seq = c(1, 5, 10), shape_seq = c(1/400, 400), n_draws = 1e4)
#'
#'
#' @export
#'
#' @importFrom stats rgamma
#' @importFrom tidyr expand_grid
#' @importFrom dplyr slice
#'
#'

generate_gaipv_prior = function(n_cats_seq,
                                shape_seq,
                                n_draws = 1e4,
                                seed = sample(.Machine$integer.max, 1)
) {
  set.seed(seed);

  n_cats_seq = sort(n_cats_seq)
  max_n_cats = max(n_cats_seq)

  all_results =
    expand_grid(n_cats = n_cats_seq,
                shape = shape_seq,
                exp_non_zero = NA)

  for(i in seq_len(nrow(all_results))) {

    curr_n_cats = slice(all_results, i) %>% pull(n_cats)
    curr_shape = slice(all_results, i) %>% pull(shape)

    base_alpha = matrix(rgamma(n_draws * (curr_n_cats + 1), shape = curr_shape), nrow = n_draws)

    curr_alpha = base_alpha[, 1:(curr_n_cats + 1), drop = FALSE]

    curr_increments = curr_alpha[,1:curr_n_cats, drop = FALSE] / rowSums(curr_alpha);

    # If the shape parameter is small, then some of the increments will be NaN
    # due to underflow. We drop those NaNs in the calculation of the mean
    all_results[i, "exp_non_zero"] = mean(rowSums(curr_increments >= 1 / (curr_n_cats + 1)), na.rm = TRUE)
    rm(curr_n_cats, curr_shape,
       curr_alpha, curr_increments)
  }
  all_results

}

