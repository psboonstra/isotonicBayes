% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_stan_divergences.R
\name{count_stan_divergences}
\alias{count_stan_divergences}
\title{Count divergences}
\usage{
count_stan_divergences(stan_fit)
}
\arguments{
\item{stan_fit}{an object of class \code{\link[rstan:stanfit-class]{rstan::stanfit()}}}
}
\value{
an integer giving the total number of divergences across all chains
}
\description{
This is a helper function to count the number of divergences
for an object of class \code{\link[rstan:stanfit-class]{rstan::stanfit()}}. It is used within \code{bayesian_isotonic}
but can also be used on its own
}
\examples{

set.seed(1)
fake_data <- make_grouped_data(runif(100), rbinom(100, 1, 0.5))
fake_model <- bayesian_isotonic(data_grouped = fake_data, return_as_stan_object = TRUE)
count_stan_divergences(fake_model)

}
\references{
\insertRef{boonstra2020b}{isotonicBayes}
}
