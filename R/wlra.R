
wlra <-
  function(x,
           w = NULL,
           lambda = 0,
           x0 = NULL,
           thresh = 1e-5,
           maxit = 1e3,
           rank.max = NULL)
  {
    d <- dim(x)
    n <- d[1]
    p <- d[2]
    if (is.null(w))
      w = matrix(rep(1, n * p), nrow = n)
    if (is.null(rank.max))
      rank.max <-
      min(n, p) - 1
    else
      rank.max <- min(rank.max, min(n, p) - 1)
    if (is.null(x0))
      x0 <- matrix(rep(0, n * p), nrow = n)
    xnas <- is.na(x)
    omega <- 1 * (!xnas)
    nz = n * p - sum(xnas)
    xfill <- x
    xfill[xnas] <- 0
    xfill <- omega * w * xfill + (1 - omega * w) * x0
    xhat <- x0
    x0 <- xfill
    iter <- 0
    error <- 100
    svd.xfill = svd(xfill)
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      svd.old = svd.xfill
      xhat.old <- xhat
      d = svd.xfill$d
      d = pmax(d - lambda, 0)
      J <- min(rank.max, length(d))
      xhat <-
        svd.xfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.xfill$v[, seq(J)]))
      xfill <- omega * w * x0 + (1 - omega * w) * xhat
      svd.xfill = svd(xfill)
      denom <-
        sum(w * (x - xhat.old) ^ 2, na.rm = T) + lambda * sum(svd.old$d[seq(J)])
      error <-
        abs(sum(w * (x - xhat) ^ 2, na.rm = T) + lambda * sum(d[seq(J)]) - denom) /
        denom
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    d <- svd.xfill$d[seq(J)]
    d = pmax(svd.xfill$d[seq(J)] - lambda, 0)
    J = min(sum(d > 0) + 1, J)
    svd.xfill = list(u = svd.xfill$u[, seq(J)], d = d[seq(J)], v = svd.xfill$v[, seq(J)])
    return(list(
      d = svd.xfill$d,
      u = svd.xfill$u,
      v = svd.xfill$v,
      cvg = cvg,
      iter = iter
    ))
  }

