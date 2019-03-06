mimi.cov <-
  function(y,
           x,
           var.type = c("gaussian", "binomial", "poisson"),
           lambda1,
           lambda2,
           algo = c("bcgd", "mcgd"),
           maxit = 100,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    yy <- y
    n <- nrow(y)
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
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n, p)-1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    q <- ncol(x)
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    alpha <- alpha0
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    theta <- theta0
    algo <- match.arg(algo,c("mcgd","bcgd"),several.ok=T)[1]
    if(algo=="bcgd"){
      res <- irwls.cov(
        y,
        x,
        var.type = var.type,
        lambda1 = lambda1,
        lambda2 = lambda2,
        maxit = maxit,
        alpha0 = alpha,
        theta0 = theta,
        thresh = thresh,
        trace.it = trace.it
      )
    } else{
      res <- mcgd.cov(y,
                      x,
                      var.type,
                      lambda1 = lambda1,
                      lambda2 = lambda2,
                      U = NULL,
                      maxit = maxit,
                      thresh = thresh,
                      alpha0 = alpha,
                      theta0 = theta,
                      R0 = NULL,
                      trace.it = trace.it)
    }

    return(list(
      y.imputed = res$y.imputed,
      param = res$param,
      alpha = res$alpha,
      theta = res$theta
    ))
  }


