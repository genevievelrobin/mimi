grad_alpha <- function(y, x, alpha, theta, var.type){
  #internal function to compute gradient of main effects coefficients
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  m <- sum(!is.na(y))
  Omega <- !is.na(y)
  q <- ncol(x)
  grad <- rep(0, n+p+q)
  par <- matrix(x%*% alpha, nrow=n) + theta
  YY <- y
  YY[is.na(y)] <- 0
  mat <- par
  mat[, var.type=="binomial"] <- exp(mat[, var.type=="binomial"])/(1+exp(mat[, var.type=="binomial"]))
  mat[, var.type=="poisson"] <- exp(mat[, var.type=="poisson"])
  mat[is.na(y)] <- 0
  grad[(n+p+1):(n+p+q)] <- t(x)%*%(-c(YY)+c(mat))/m
  return(grad)
}
