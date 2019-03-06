mimi.lr <-
  function(y,
           var.type = c("gaussian", "binomial", "poisson"),
           lambda1,
           algo = c("bcgd", "mcgd"),
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    yy <- y
    if (sum(var.type == "binomial") > 0) {
      for (j in 1:sum(var.type == "binomial")) {
        y <- data.frame(y)
        y[, which(var.type == "binomial")[j]] <-
          as.factor(y[, which(var.type == "binomial")[j]])
      }
    }
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n,p)-1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    theta <- theta0
    algo <- match.arg(algo,c("mcgd","bcgd"),several.ok=T)[1]
    if(algo == "bcgd"){
      res <-
        irwls.lr(
          y,
          var.type = var.type,
          lambda1 = lambda1,
          maxit = maxit,
          theta0 = theta,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
    } else {
      res <- mcgd.lr(y,
                     var.type,
                     lambda1 = lambda1,
                     maxit = maxit,
                     thresh = thresh,
                     theta0 = theta0,
                     trace.it = trace.it)

    }
    return(list(y.imputed = res$y.imputed, theta = res$theta))
  }


