#' @title Categorize continuous data
#'
#' @description This is a helper function to turn a vector of continuous predictors x
#' and corresponding vector of binary outcomes y into a data.frame of a form
#' expected by the function [bayesian_isotonic()]
#'
#' @param x A vector of continuous predictors
#'
#' @param y A vector of Bernoulli outcomes with length equal to length(x)
#'
#' @param breaks A vector of values used to categorize x
#'
#'
#' @references
#' \insertRef{boonstra2020b}{isotonicBayes}
#'
#' @examples
#'
#' @importFrom dplyr %>% group_by summarize arrange mutate
#' @importFrom stats quantile
#' @importFrom tibble tibble
#' @importFrom Rdpack reprompt
#'

make_grouped_data <- function(x, y, breaks = NULL) {

  if(is.null(breaks)) {
    message("Categorizing 'x' into five groups. Consider choosing your own value.")
    breaks <-
      quantile(x, probs = seq(0, 1, length = 6)) %>%
      as.numeric() %>%
      replace(which.min(.), -Inf) %>%
      replace(which.max(.), Inf)
  }

  tibble(x = cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE),
         y = y) %>%
    group_by(x) %>%
    summarize(n = length(x),
              y = sum(y)) %>%
    arrange(x) %>%
    mutate(x_cat = 1:n())
}
