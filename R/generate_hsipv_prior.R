#' @title Generate draws from HSIPV prior
#'
#' @description Generate draws from the HSIPV prior and then calculate the
#' effective number of parameters according to this prior.
#'
#' @param n_cats_seq A vector of one or more positive integers giving the
#' desired values of K over which to evaluate
#'
#' @param scale_seq A vector of one or more positive numbers giving the
#' set of scale parameters c over which to evaluate
#'
#'
#' @param local_dof a positive number giving the degrees of freedom for the
#' t-distribution of the local shrinkage parameter.
#' local_dof = global_dof = 1 corresponds to the horseshoe
#'
#' @param global_dof a positive number giving the degrees of freedom for the
#' t-distribution of the global shrinkage parameter.
#' local_dof = global_dof = 1 corresponds to the horseshoe
#'
#'
#' @param slab_precision a positive number giving the precision of the slab
#' parameter for the regularized horseshoe. Defaults to 1, which is the only
#' value considered in the manuscript
#'
#'
#' @param n_draws a positive integer giving the number of random draws from the
#' t-distributions with which to estimate expectations
#'
#' @param seed a positive integer giving the random seed, provided for reproducibility.
#'
#' @return a tibble with number of rows equal to the product of the lengths of
#' `n_cats_seq` and `scale_seq` and three columns: the current number of categories,
#' the current scale, and the calculated effective number of prior parameters
#'
#' @references
#'
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' \insertRef{piironen2017hyperprior}{isotonicBayes}
#'
#' @examples
#' generate_hsipv_prior(n_cats_seq = c(1, 5, 10), scale_seq = c(1/400, 400), n_draws = 1e4)
#'
#'
#' @export
#'
#' @importFrom stats rt rnorm
#' @importFrom matrixStats rowCumsums
#' @importFrom tidyr expand_grid
#' @importFrom dplyr slice
#'
#'

generate_hsipv_prior = function(n_cats_seq,
                                scale_seq,
                                local_dof = 1,
                                global_dof = 1,
                                slab_precision = 1,
                                n_draws,
                                seed = sample(.Machine$integer.max, 1)
) {
  set.seed(seed);

  n_cats_seq = sort(n_cats_seq)
  max_n_cats = max(n_cats_seq)

  lambda_tau_sq =
    matrix(
      rt(n_draws * (max_n_cats + 1), df = local_dof)^2 *
        rt(n_draws, df = global_dof)^2,
      nrow = n_draws)

  base_alpha = matrix(abs(rnorm(n_draws * (max_n_cats + 1))), nrow = n_draws);
  all_results =
    expand_grid( n_cats = n_cats_seq, scale = scale_seq,exp_non_zero = NA)

  for(i in seq_len(nrow(all_results))) {
    curr_n_cats = slice(all_results, i) %>% pull(n_cats)
    curr_scale = slice(all_results, i) %>% pull(scale)
    curr_base_alpha = base_alpha[, 1:(curr_n_cats + 1), drop = FALSE]

    foo = lambda_tau_sq[, 1:(curr_n_cats + 1), drop = FALSE] * curr_scale^2
    curr_prior_scales = sqrt(1 / (slab_precision + 1 / foo));
    curr_alpha = curr_prior_scales * curr_base_alpha
    curr_increments = curr_alpha[,1:curr_n_cats, drop = FALSE] / rowSums(curr_alpha);

    all_results[i, 3] = mean(rowSums(curr_increments >= 1 / (curr_n_cats + 1)))
    rm(curr_n_cats, curr_scale, curr_base_alpha, foo,
       curr_prior_scales, curr_alpha, curr_increments)
  }
 all_results

}

