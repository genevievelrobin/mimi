quad_approx <- function(y, param, var.type) {
  n <- nrow(y)
  p <- length(var.type)
  yt <- rep(0, n)
  vt <- rep(0, n)
  for (j in 1:p) {
    w <-
      lapply(1:n, function(i)
        wght(
          y = y[i, j],
          param = param[i, j],
          var.type = var.type[j]
        ))
    ytilde <- sapply(1:n, function(i)
      w[[i]]$ytilde)
    vtilde2 <- sapply(1:n, function(i)
      w[[i]]$vtilde2)
    yt <- cbind(yt, as.matrix(ytilde))
    vt <- cbind(vt, as.matrix(vtilde2))
  }
  return(list(ytilde = yt[, 2:ncol(yt)], vtilde2 = vt[, 2:ncol(vt)]))
}
