#' Turning variable into vector if scalar.
#' 
#' @param scalar The variable to be turned into a vector if not of length reps.
#' @param reps The length that the variable should have.
#' 
#' @return A vector of length reps, with scalar repeated until the appropriate
#' length is achieved.
ScalarToVector <- function(scalar, reps) {
  # Functions that repeats a scalar multiple times.
  if (length(scalar) < reps) {
    scalar <- rep(scalar, reps)
  }
  return(scalar)
}