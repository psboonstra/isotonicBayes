test_that("horseshooe_good", {
  expect_silent(bayesian_isotonic(data_grouped = tibble(y = c(0,1,2), n = c(2,2,2)),
                                  prior_type = "horseshoe",
                                  stan_args = list(global_dof_stan = 1,
                                                   local_dof_stan = 1,
                                                   alpha_scale_stan = 1),
                                  stan_seed = 1))
})

test_that("gamma_good", {
  expect_silent(bayesian_isotonic(data_grouped = tibble(y = c(0,1,2), n = c(2,2,2)),
                                  prior_type = "gamma",
                                  stan_args = list(alpha_shape_stan = 0.1,
                                                   tiny_positive_stan = 1e-5),
                                  stan_seed = 1))
})

