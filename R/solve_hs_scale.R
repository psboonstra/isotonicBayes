#' solve_hs_scale
#'
#' @param target_mean
#'
#'
#' @references
#' \insertRef{boonstra2020}{isotonicBayes}
#'
#' @examples
#'
#' @importFrom stats rt rnorm
#' @importFrom Rdpack reprompt
#'



solve_for_hs_scale = function(target_mean,
                              local_dof = 1,
                              global_dof = 1,
                              slab_precision = 1,
                              n,
                              sigma = 2,
                              tol = .Machine$double.eps,
                              max_iter = 100,
                              n_sim = 1e6
) {

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
  abs_normal_draws = abs(rnorm(n_sim));
  stopifnot(target_mean > 0 && target_mean < 1);#Ensure proper bounds
  log_scale = diff_target = numeric(max_iter);
  log_scale[1] = log(target_mean/(1 - target_mean)*sigma/sqrt(n));
  prior_scales = 1 / sqrt(slab_precision + 1/(exp(2*log_scale[1]) * lambda^2));
  kappa = 1/(1+n*prior_scales/sigma^2);
  diff_target[1] = mean(1-kappa) - target_mean;
  log_scale[2] = 0.5 + log_scale[1];
  prior_scales = 1 / sqrt(slab_precision + 1/(exp(2*log_scale[2]) * lambda^2));
  kappa = 1/(1+n*prior_scales/sigma^2);
  diff_target[2] = mean(1-kappa) - target_mean;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_scale[i] =
      log_scale[i-1] -
      diff_target[i-1]*(log_scale[i-1]-log_scale[i-2])/(diff_target[i-1]-diff_target[i-2]);
    prior_scales =
      1 / sqrt(slab_precision + 1/(exp(2*log_scale[i]) * lambda^2));
    kappa = 1/(1+n*prior_scales/sigma^2);
    diff_target[i] = mean(1-kappa) - target_mean;
    if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
  }

  list(scale = exp(log_scale[i]),
       achieved_mean = mean(1-kappa),
       target_mean = target_mean,
       diff_from_target = abs(diff_target[i]),
       iter = i);
}

