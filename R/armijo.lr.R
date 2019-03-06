armijo.lr <-
  function(y0,
           theta,
           theta.tmp,
           alpha.mat=0,
           w,
           z,
           lambda1,
           var.type,
           b = 0.5,
           thresh = 1e-5,
           nu=0.1,
           zeta=0.1,
           th=0.1,
           step=1) {
    d <- dim(y0)
    n <- d[1]
    p <- d[2]
    param.tmp <- theta.tmp+alpha.mat
    omega <- !is.na(y0)
    gaus.tmp <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois.tmp <-
      sum(-(y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
            exp(param.tmp[, var.type == "poisson"]), na.rm = T)
    binom.tmp <-
      sum(-(y0[, var.type == "binomial"] * param.tmp[, var.type == "binomial"]) +
            log(1 + exp(param.tmp[, var.type == "binomial"])), na.rm = T)
    d.tmp <- svd(theta.tmp)$d
    flag <- TRUE
    direction <- theta
    d0 <- svd(theta.tmp + direction)$d
    while (flag) {
      theta <- theta.tmp + step * direction
      param <- theta+alpha.mat
      gaus <-
        (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                        2, na.rm = T)
      pois <-
        sum(-(y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
              exp(param[, var.type == "poisson"]), na.rm = T)
      binom <-
        sum(-(y0[, var.type == "binomial"] * param[, var.type == "binomial"]) +
              log(1 + exp(param[, var.type == "binomial"])), na.rm = T)
      d <- svd(theta)$d
      diff <- pois-pois.tmp + gaus-gaus.tmp + binom- binom.tmp + lambda1 * (sum(d) - sum(d.tmp))
      if (diff <= step * zeta * (-2 * sum(omega*w * z * direction, na.rm = T) +
                                 th * sum((w*omega + (1-omega)*nu) * direction ^ 2) + lambda1 * (sum(d0) - sum(d.tmp)))) {
        flag <- FALSE
      } else if (step<=1e-30){
        flag <- FALSE
      } else{
        step <- step / 2
      }
    }
    obj <- pois + gaus + binom + lambda1 * d
    return(list(
      theta = theta,
      objective = obj,
      step = step
    ))
  }
