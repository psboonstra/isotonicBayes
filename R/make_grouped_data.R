
make_grouped_data <- function(x, y, breaks = NULL) {

  if(is.null(breaks)) {
    message("Categorizing 'x' into five groups. Consider choosing your own value.")
    breaks <-
      quantile(x, probs = seq(0, 1, length = 6)) %>%
      as.numeric() %>%
      replace(which.min(.), -Inf) %>%
      replace(which.max(.), Inf)
  }

  tibble(x = cut(x, breaks = breaks, right = F),
         y = y) %>%
    group_by(x) %>%
    summarize(n = length(x),
              y = sum(y)) %>%
    arrange(x) %>%
    mutate(x_cat = 1:n())
}
