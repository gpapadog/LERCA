ScalarToVector <- function(scalar, reps) {
  # Functions that repeats a scalar multiple times.
  if (length(scalar) < reps) {
    scalar <- rep(scalar, reps)
  }
  return(scalar)
}