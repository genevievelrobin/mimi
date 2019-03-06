armijo.alpha <-
  function(y0,
           x,
           alpha,
           theta,
           alpha.tmp,
           theta.tmp,
           z,
           w,
           b = 0.5,
           lambda2,
           var.type,
           thresh = 1e-5,
           nu=0.1,
           zeta=0.9,
           th=0.9,
           step=1) {
    d <- dim(y0)
    n <- d[1]
    p <- d[2]
    q <- ncol(x)
    omega <- !is.na(y0)
    direction <- alpha
    direction.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    alpha.tmp.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha.tmp, nrow = n)
    param.tmp <- alpha.tmp.mat + theta.tmp
    gaus.tmp <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois.tmp <-
      sum((-(y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
             exp(param.tmp[, var.type == "poisson"])), na.rm = T)
    binom.tmp <-
      sum((-(y0[, var.type == "binomial"] * param.tmp[, var.type == "binomial"]) +
             log(1 + exp(param.tmp[, var.type == "binomial"]))), na.rm = T)
    flag <- TRUE
    a0 <- sum(abs(alpha.tmp + direction))
    while (flag) {
      alpha <- alpha.tmp + step * direction
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      param <- alpha.mat + theta
      gaus <-
        (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum((-(y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
               exp(param[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum((-(y0[, var.type == "binomial"] * param[, var.type == "binomial"]) +
               log(1 + exp(param[, var.type == "binomial"]))), na.rm = T)
      diff <- gaus - gaus.tmp + pois - pois.tmp + binom - binom.tmp + lambda2 * (sum(abs(alpha)) - sum(abs(alpha.tmp)))
      if (diff <= step * zeta * (-2 * sum(omega*w * z * alpha.tmp.mat) + th *
                                 sum((omega*w + (1-omega)*nu) * alpha.tmp.mat ^ 2) + lambda2 * (a0 - sum(abs(alpha.tmp))))) {
        flag <- FALSE
      } else if (step <=1e-30){
        flag <- FALSE
      } else
        step <- b * step
    }
    return(list(alpha = alpha,
                step=step))
  }
