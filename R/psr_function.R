#' @export
psr <- function(samples) {
  R <- nrow(samples)
  if (any(is.na(samples))) {
    warning('Some values are NA. na.rm is set to TRUE.') 
  }
  B <- R * var(colMeans(samples, na.rm = TRUE))
  W <- mean(apply(samples, 2, var, na.rm = TRUE))
  value <- sqrt((B + (R - 1) * W) / (R * W))
  return(value)
}
