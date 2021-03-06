% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_isotonic.R
\name{bayesian_isotonic}
\alias{bayesian_isotonic}
\title{Bayesian isotonic regression}
\usage{
bayesian_isotonic(
  data_grouped = NULL,
  prior_type = "horseshoe",
  stan_args = list(local_dof_stan = 1, global_dof_stan = 1, alpha_scale_stan = 1),
  sample_from_prior_only = FALSE,
  conf_level = 0.5,
  conf_level_direction = "both",
  sig_threshold = c(0.005, 0.01, 0.05),
  verbose = FALSE,
  mc_warmup = 2500,
  mc_samps = 5000,
  mc_chains = 1,
  mc_thin = 1,
  mc_stepsize = 0.1,
  mc_adapt_delta = 0.99,
  mc_max_treedepth = 15,
  return_as_stan_object = FALSE,
  stan_seed = sample.int(.Machine$integer.max, 1)
)
}
\arguments{
\item{data_grouped}{a tibble or data.frame with named columns y and n indicating,
respectively, the number of events and number of bernoulli trials. Each
row corresponds to a different category of the predictor. \strong{The rows are assumed
to be ordered such that the first row corresponds to the first
category of the predictor, the second row the second category, and so forth}}

\item{prior_type}{a character indicating the type of prior to use. Currently
can be only "horseshoe" or "gamma"}

\item{stan_args}{a named list of arguments corresponding to the selected prior.
If prior_type="horseshoe", then stan_args must contain the following named
components: local_dof_stan (an integer giving the degrees of freedom for
the local shrinkage parameter), global_dof_stan (same for the global
shrinkage parameter), and alpha_scale_stan (the scale tuning parameter).
If prior_type="gamma", then stan_args must contain the following named
components: alpha_shape_stan (the shape tuning parameter) and tiny_positive_stan
(a small number specifying the lower truncation of the gamma distribution)}

\item{sample_from_prior_only}{a logical that offers an easy way to sample from the
prior. If TRUE, the data provided via data_grouped are ignored}

\item{conf_level}{a number between 0 and 1 specifying the posterior quantiles
to be calculated. The default value is 0.50, meaning that the function
will by default return the 25th and 75th percentiles of the posterior
distribution of each parameter}

\item{conf_level_direction}{a character that must be "both", "lower", or "upper"}

\item{sig_threshold}{a vector of proportions. For each proportion, the posterior
probability that each increment exceeds that proportion will be calculated
and returned}

\item{verbose}{a logical indicating whether to return all posterior draws
of all parameters}

\item{mc_warmup}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{warmup}}

\item{mc_samps}{the sum of \code{mc_warmup} and \code{mc_samps} is passed to
\code{\link[rstan:stan]{stan()}} as the value of \code{iter}}

\item{mc_chains}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{chains}}

\item{mc_thin}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{thin}}

\item{mc_stepsize}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{stepsize}
(inside \code{control})}

\item{mc_adapt_delta}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{adapt_delta}
(inside \code{control})}

\item{mc_max_treedepth}{passed to \code{\link[rstan:stan]{stan()}} as the value of \code{max_treedepth}
(inside \code{control})}

\item{return_as_stan_object}{a logical. Should the function return an object
of class \code{stanfit}? Defaults to \code{FALSE}}

\item{stan_seed}{A positive integer to seed}
}
\description{
bayesian_isotonic
}
\examples{

set.seed(1)
fake_data <- make_grouped_data(runif(100), rbinom(100, 1, 0.5))
bayesian_isotonic(data_grouped = fake_data)

}
\references{
\insertRef{boonstra2020b}{isotonicBayes}
}
