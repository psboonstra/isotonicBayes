#' @title Solve for m_eff / (K - 1)
#'
#' @description Takes as input a scale parameter, interpreted as c from in Equation (2.9),
#' and returns m_eff / (K - 1),  where m_eff is also defined in Equation (2.9) and K is
#' the number of categories.
#'
#' @param scale a positive number corresponding to the scale parameter
#'
#' @param local_dof a positive number giving the degrees of freedom for the
#' t-distribution of the local shrinkage parameter.
#' local_dof = global_dof = 1 corresponds to the horseshoe
#'
#' @param global_dof a positive number giving the degrees of freedom for the
#' t-distribution of the global shrinkage parameter.
#' local_dof = global_dof = 1 corresponds to the horseshoe
#'
#' @param slab_precision a positive number giving the precision of the slab
#' component of the regularized horseshoe.
#'
#' @param n a positive integer giving the sample size
#'
#' @param sigma the presumed standard deviation parameter for the measurement error.
#' If the outcome is binary, then sigma=2 corresponds to the most conservative
#' choice in the sense of yielding more shrinkage than desired \insertCite{piironen2017hyperprior}{isotonicBayes}
#'
#' @param tol a small number giving the acceptable numerical tolerance
#' for determining equality
#'
#' @param max_iter a positive integer giving the maximum number of iterations
#' to proceed without convergence before giving up
#'
#' @param n_sim a positive integer giving the number of random draws from the
#' t-distributions with which to estimate expectations
#'
#' @param seed a positive integer giving the random seed, provided for reproducibility.
#'
#' @references
#'
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' \insertRef{piironen2017hyperprior}{isotonicBayes}
#'
#'
#' @importFrom stats rt rnorm var
#' @importFrom Rdpack reprompt
#'



solve_for_hs_propeff = function(scale,
                                local_dof = 1,
                                global_dof = 1,
                                slab_precision = 1,
                                n,
                                sigma = 2,
                                tol = .Machine$double.eps,
                                max_iter = 100,
                                n_sim = 1e6,
                                seed = sample(.Machine$integer.max, 1)
) {
  set.seed(seed);
  stopifnot(slab_precision >= 0 && n > 0 && sigma > 0);#Ensure proper bounds
  stopifnot(local_dof >= 0 && global_dof >= 0);#Ensure proper bounds
  stopifnot(scale > 0);#Ensure proper bounds
  do_local = (local_dof > 0);
  do_global = (global_dof > 0);
  if(do_local) {
    lambda = rt(n_sim,df = local_dof);
  } else {
    lambda = rep(1, n_sim);
  }
  if(do_global) {
    lambda = lambda * rt(n_sim,df = global_dof);
  }
  prior_scales_sq = 1 / (slab_precision + 1 / (scale^2 * lambda^2));
  kappa = 1 / (1 + n * prior_scales_sq / sigma^2);

  list(mean_prop = mean(1 - kappa),
       var_prop = var(1 - kappa))
}

