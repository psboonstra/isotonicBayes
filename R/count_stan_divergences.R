

count_stan_divergences = function(stan_fit) {
  foo = get_sampler_params(stan_fit, inc_warmup = FALSE);
  n_draws = lapply(foo, nrow)[[1]];
  sum(unlist(lapply(foo,"[",i = 1:n_draws, j = "divergent__")));
}
