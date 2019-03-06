log_factorial <- function(x) {
  m <- apply(x, c(1, 2), function(t)
    if (is.na(t)) {
      NA
    } else if (t == 0) {
      1
    } else
      sum(log(1:t)))
  return(m)
}
