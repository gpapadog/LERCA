get_stats <- function(x, probs) {
  x <- as.numeric(x)
  quants <- quantile(x, probs = probs)
  names(quants) <- NULL
  return(c(mean = mean(x), LB = quants[1], UB = quants[2]))
}