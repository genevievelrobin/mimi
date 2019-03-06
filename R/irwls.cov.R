irwls.cov <-
  function(y,
           x,
           var.type,
           lambda1,
           lambda2,
           maxit = 100,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL,
           nu = 1e-2) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if (is.null(max.rank))
      max.rank <- min(n, p) - 1
    else
      max.rank <- min(max.rank, min(n, p) - 1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    y0 <- y
    q <- ncol(x)
    x <- as.matrix(x)
    x <- matrix(as.numeric(x), nrow = nrow(x))
    omega <- !is.na(y)
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    alpha <- alpha0
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    theta <- theta0
    param <- alpha.mat + theta
    gaus <- (1 / 2) * sum(y0[, var.type == "gaussian"] ^ 2, na.rm = T)
    pois <- n * p
    binom <- binom <- -log(2) * n * p
    objective <- gaus + pois + binom
    y0 <- y
    error <- 1
    iter <- 0
    stepalpha <- 1
    steptheta <- 1
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      alpha.tmp <- alpha
      theta.tmp <- theta
      yv <- quad_approx(y0, param, var.type)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega == 0] <- 0
      alpha <-
        glmnet::glmnet(
          x,
          c((ytilde - omega * theta.tmp) * vtilde2 / (omega * vtilde2 + (1 -
                                                                           omega) * nu)),
          family = "gaussian",
          lambda = lambda2,
          intercept = FALSE,
          weights = c(omega * vtilde2 + (1 - omega) * nu)
        )
      alpha <- as.numeric(alpha$beta)
      resalpha <-
        armijo.alpha(
          y,
          x,
          alpha,
          theta.tmp,
          alpha.tmp,
          theta.tmp,
          z = ytilde - omega * param,
          w = vtilde2,
          lambda2 = lambda2,
          var.type = var.type,
          nu = nu,
          zeta = 0.1,
          th = 0,
          step = stepalpha
        )
      alpha <- resalpha$alpha
      stepalpha <- resalpha$step
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      alpha.tmp <- alpha
      param <- alpha.mat + theta.tmp
      yv <- quad_approx(y0, param, var.type)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega == 0] <- 0
      w <- omega * vtilde2 + (1 - omega) * nu
      lambda1w <- lambda1 / max(w)
      w <- w / max(w)
      svd_theta <-
        wlra(vtilde2 / (omega * vtilde2 + (1 - omega) * nu) * (ytilde - omega * alpha.mat),
             w,
             lambda1w)
      u <- svd_theta$u
      d <- svd_theta$d
      v <- svd_theta$v
      if (is.null(dim(u))) {
        theta <- d * u %*% t(v)
      } else {
        theta <- u %*% diag(d) %*% t(v)
      }
      restheta <-
        armijo.lr(
          y0 = y,
          theta = theta,
          theta.tmp = theta.tmp,
          alpha.mat = alpha.mat,
          w = vtilde2,
          z = ytilde - omega * param,
          lambda1 = lambda1,
          var.type = var.type,
          b = 0.5,
          nu = nu,
          zeta = 0.1,
          th = 0,
          step = steptheta
        )
      theta <- restheta$theta
      steptheta <- restheta$step
      param <- alpha.mat + theta
      gaus <-
        (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum((-(y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
               exp(param[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum((-(y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
               log(1 + exp(param[, var.type == "binomial"]))), na.rm = T)
      d <- svd(theta)$d
      objective <- c(objective, min(
        .Machine$double.xmax,
        (pois + gaus + binom + lambda1 * sum(d) + lambda2 * sum(abs(alpha)))
      ))
      if (iter == 1) {
        error <- 1
      } else if ((stepalpha <= 1e-30) && (steptheta <= 1e-30)) {
        error <- 0
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <-
            abs(objective[iter] - objective[iter - 1]) / abs(objective[iter - 1])

      }
      if (trace.it) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    estim <- param
    estim[, var.type == "poisson"] <-
      matrix(stats::rpois(n * sum(var.type == "poisson"), lambda = c(exp(estim[, var.type == "poisson"]))), nrow =
               n)
    estim[, var.type == "binomial"] <-
      round(exp(estim[, var.type == "binomial"]) / (1 + exp(estim[, var.type == "binomial"])))
    y.imputed <- y0
    y.imputed[omega == 0] <- estim[omega == 0]
    return(list(
      y.imputed = y.imputed,
      param = alpha.mat + theta,
      alpha = alpha,
      theta = theta
    ))
  }
