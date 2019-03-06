irwls.lr <-
  function(y,
           var.type,
           lambda1,
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-6,
           trace.it = F,
           max.rank = NULL,
           nu=0.1) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n,p)-1)
    y <- matrix(as.numeric(as.matrix(y)), nrow = n)
    y0 <- y
    omega <- !is.na(y)
    if (is.null(theta0)) {
      theta0 <- matrix(rep(0, n * p), nrow = n)
    }
    theta <- theta0
    gaus <- (1 / 2) * sum(y0[, var.type == "gaussian"]^2, na.rm = T)
    pois <- n*p
    binom <--log(2)*n*p
    objective <- gaus+pois+binom

    error <- 1
    iter <- 0
    step <- 1
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      theta.tmp <- theta
      yv <- quad_approx(y0, theta.tmp, var.type)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega == 0] <- 0
      z <- ytilde * vtilde2 / (vtilde2 + nu)
      w <- vtilde2*omega + (1-omega)*nu
      lambda1w <- lambda1 / max(w)
      w <- w / max(w)
      svd_theta <-
        wlra(
          x = z,
          w = w,
          lambda = lambda1w,
          x0 = NULL,
          thresh = 0.1 * thresh,
          rank.max = max.rank
        )
      u <- svd_theta$u
      d <- svd_theta$d
      v <- svd_theta$v
      if (is.null(dim(u))) {
        theta <- d * u %*% t(v)
      } else {
        theta <- u %*% diag(d) %*% t(v)
      }
      direction_theta <- theta
      res <-
        armijo.lr(
          y0=y0,
          theta=theta,
          theta.tmp=theta.tmp,
          alpha.mat=0,
          w = vtilde2,
          z = ytilde - omega*theta.tmp,
          lambda1=lambda1,
          var.type=var.type,
          b = 0.5,
          thresh=thresh,
          nu=nu,
          zeta=0.1,
          th=0.1,
          step=step
        )
      step <- res$step
      theta <- res$theta
      objective <-
        c(objective, min(.Machine$double.xmax, res$objective))
      if (iter == 1) {
        error <- 1
      } else if(step<=1e-30){
        error <- 0
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <- abs(objective[iter]-objective[iter-1])/abs(objective[iter-1])
      }
      if(all(var.type=="gaussian")){
        error <- 0
      }
      if (trace.it) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    estim <- theta
    estim[, var.type == "poisson"] <- exp(estim[, var.type == "poisson"])
    estim[, var.type == "binomial"] <- round(exp(estim[, var.type == "binomial"]) / (1 + exp(estim[, var.type == "binomial"])))
    y.imputed <- y0
    y.imputed[omega==0] <- estim[omega==0]
    return(list(
      theta = theta,
      y.imputed = y.imputed
    ))
  }
