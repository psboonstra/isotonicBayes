#' @title Bayesian isotonic regression
#'
#' @description  bayesian_isotonic
#'
#' @param data_grouped a tibble or data.frame with named columns y and n indicating,
#' respectively, the number of events and number of bernoulli trials. Each
#' row corresponds to a different category of the predictor. **The rows are assumed
#' to be ordered such that the first row corresponds to the first
#' category of the predictor, the second row the second category, and so forth**
#'
#' @param prior_type a character indicating the type of prior to use. Currently
#' can be only "horseshoe" or "gamma"
#'
#' @param stan_args a named list of arguments corresponding to the selected prior.
#' If prior_type="horseshoe", then stan_args must contain the following named
#' components: local_dof_stan (an integer giving the degrees of freedom for
#' the local shrinkage parameter), global_dof_stan (same for the global
#' shrinkage parameter), and alpha_scale_stan (the scale tuning parameter).
#' If prior_type="gamma", then stan_args must contain the following named
#' components: alpha_shape_stan (the shape tuning parameter) and tiny_positive_stan
#' (a small number specifying the lower truncation of the gamma distribution)
#'
#' @param sample_from_prior_only a logical that offers an easy way to sample from the
#' prior. If TRUE, the data provided via data_grouped are ignored
#'
#' @param conf_level a number between 0 and 1 specifying the posterior quantiles
#' to be calculated. The default value is 0.50, meaning that the function
#' will by default return the 25th and 75th percentiles of the posterior
#' distribution of each parameter
#'
#' @param conf_level_direction a character that must be "both", "lower", or "upper"
#'
#' @param sig_threshold a vector of proportions. For each proportion, the posterior
#' probability that each increment exceeds that proportion will be calculated
#' and returned
#'
#' @param verbose a logical indicating whether to return all posterior draws
#' of all parameters
#'
#' @param mc_warmup passed to \code{\link[rstan:stan]{stan()}} as the value of `warmup`
#'
#' @param mc_samps the sum of `mc_warmup` and `mc_samps` is passed to
#' \code{\link[rstan:stan]{stan()}} as the value of `iter`
#'
#' @param mc_chains passed to \code{\link[rstan:stan]{stan()}} as the value of `chains`
#'
#' @param mc_thin passed to \code{\link[rstan:stan]{stan()}} as the value of `thin`
#'
#' @param mc_stepsize passed to \code{\link[rstan:stan]{stan()}} as the value of `stepsize`
#' (inside `control`)
#'
#' @param mc_adapt_delta passed to \code{\link[rstan:stan]{stan()}} as the value of `adapt_delta`
#' (inside `control`)
#'
#' @param mc_max_treedepth passed to \code{\link[rstan:stan]{stan()}} as the value of `max_treedepth`
#' (inside `control`)
#'
#' @param return_as_stan_object a logical. Should the function return an object
#' of class `stanfit`? Defaults to `FALSE`
#'
#'
#' @param stan_seed A positive integer to seed
#'
#' @references
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' @examples
#'
#' set.seed(1)
#' fake_data <- make_grouped_data(runif(100), rbinom(100, 1, 0.5))
#' bayesian_isotonic(data_grouped = fake_data)
#'
#' @export
#'
#' @importFrom dplyr %>% bind_cols mutate pull
#' @importFrom tibble as_tibble
#' @importFrom rlang := sym
#' @importFrom rstan stan get_elapsed_time sampling
#' @importFrom Rdpack reprompt
#'

bayesian_isotonic = function(data_grouped = NULL,
                             prior_type = "horseshoe",
                             stan_args = list(
                               local_dof_stan = 1,
                               global_dof_stan = 1,
                               alpha_scale_stan = 1),
                             sample_from_prior_only = FALSE,
                             conf_level = 0.50,
                             conf_level_direction = "both",
                             sig_threshold = c(0.005, 0.01, 0.05),
                             verbose = FALSE,
                             mc_warmup = 2.5e3,
                             mc_samps = 5e3,
                             mc_chains = 1,
                             mc_thin = 1,
                             mc_stepsize = 0.1,
                             mc_adapt_delta = 0.99,
                             mc_max_treedepth = 15,
                             return_as_stan_object = FALSE,
                             stan_seed = sample.int(.Machine$integer.max, 1)) {

  stopifnot(isTRUE("data.frame" %in% class(data_grouped)))
  stopifnot(isTRUE(all(c("y","n") %in% colnames(data_grouped))))
  stopifnot(all(pull(data_grouped,y) >= 0) &&
              all(pull(data_grouped,y) <= pull(data_grouped,n)));
  stopifnot(conf_level >= 0 && conf_level <= 1);
  stopifnot(prior_type %in% c("horseshoe", "gamma",
                              "horseshoe2", "gamma2",
                              "horseshoe3"))

  if(prior_type == "horseshoe" && !setequal(names(stan_args), c("local_dof_stan", "global_dof_stan", "alpha_scale_stan")))
    stop("The horseshoe prior expects that stan_args must contain all and only named elements local_dof_stan, global_dof_stan, and alpha_scale_stan")

  if(prior_type == "gamma" && !setequal(names(stan_args), c("alpha_shape_stan", "tiny_positive_stan")))
    stop("The gamma prior expects that stan_args must contain all and only named elements alpha_shape_stan and tiny_positive_stan")



  curr_fit <- sampling(object = stanmodels[[paste0("iso_",prior_type)]],
                       data = c(list(n_groups_stan = nrow(data_grouped),
                                     n_per_group_stan = as.array(pull(data_grouped,n)),
                                     y_stan = as.array(pull(data_grouped,y)),
                                     only_prior_stan = as.integer(sample_from_prior_only)),
                                stan_args),
                       warmup = mc_warmup,
                       iter = mc_samps + mc_warmup,
                       chains = mc_chains,
                       thin = mc_thin,
                       control = list(adapt_delta = mc_adapt_delta,
                                      stepsize = mc_stepsize,
                                      max_treedepth = mc_max_treedepth),
                       seed = stan_seed,
                       verbose = F,
                       refresh = 0)

  number_divergences = count_stan_divergences(curr_fit);
  max_rhat = max(rstan::summary(curr_fit)$summary[,"Rhat"], na.rm = TRUE)

  if(!return_as_stan_object) {
    chain_run_times_secs = rowSums(get_elapsed_time(curr_fit));
    total_run_time_secs = max(chain_run_times_secs);
    foo = rstan::extract(curr_fit);

    xi_number_nan = colSums(is.na(foo$xi))
    alpha_number_nan = colSums(is.na(foo$alpha))
    mean_prob = colMeans(foo$xi, na.rm = T);

    if(conf_level_direction == "both") {
      quantile_probs <-
        apply(foo$xi,
              2,
              quantile,
              p = 1/2 + c(-conf_level, 0, conf_level)/2, na.rm = T);
    } else if(conf_level_direction == "lower") {
      quantile_probs <-
        rbind(apply(foo$xi,
                    2,
                    quantile,
                    p = c(1 - conf_level, 1/2), na.rm = T),
              "100%" = 1);
    } else {
      quantile_probs <-
        rbind("0%" = 0,
              apply(foo$xi,
                    2,
                    quantile,
                    p = c(1/2, conf_level), na.rm = T));
    }
    quantile_probs = t(quantile_probs);
    colnames(quantile_probs) = c("model_lower_ci_prob", "model_median_prob", "model_upper_ci_prob");
  } else {
    foo = rstan::extract(curr_fit);
    xi_number_nan = colSums(is.na(foo$xi))
    alpha_number_nan = colSums(is.na(foo$alpha))
  }



  if(number_divergences > 0) {
    warning(paste0("there were ", number_divergences, " divergent transitions"));
  }
  if(any(xi_number_nan > 0) | any(alpha_number_nan > 0))  {
    warning(paste0("there were ", max(c(xi_number_nan,alpha_number_nan)), " draws in which one or more elements of xi were NaN"));
  }

  if(return_as_stan_object) {
    curr_fit;
  } else {

    data_grouped =
      data_grouped %>%
      mutate(emp_mean_prob = y/n) %>%
      bind_cols(model_mean_prob = mean_prob,
                as_tibble(quantile_probs));

    draws_delta = t(apply(cbind(0,0,foo$xi),1,diff))[,-1,drop = F];
    for(i in 1:length(sig_threshold)) {
      data_grouped =
        bind_cols(data_grouped,
                  !!sym(paste0("prob_delta_gt_",sig_threshold[i])) := colMeans(draws_delta > sig_threshold[i], na.rm = T));
    }

    if(verbose) {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level,
             prior_type = prior_type),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             number_divergences = number_divergences,
             max_rhat = max_rhat,
             xi_number_nan = xi_number_nan,
             alpha_number_nan = alpha_number_nan,
             any_nan = max(c(xi_number_nan, alpha_number_nan) > 0),
             all_draws = foo,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));
    } else {
      c(list(data_grouped = data_grouped,
             conf_level = conf_level,
             prior_type = prior_type),
        stan_args,
        list(sample_from_prior_only = sample_from_prior_only,
             number_divergences = number_divergences,
             max_rhat = max_rhat,
             xi_number_nan = xi_number_nan,
             alpha_number_nan = alpha_number_nan,
             any_nan = max(c(xi_number_nan, alpha_number_nan) > 0),
             all_draws = NA,
             chain_run_times_secs = chain_run_times_secs,
             total_run_time_secs = total_run_time_secs));

    }
  }
}
