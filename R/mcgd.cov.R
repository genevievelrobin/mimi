mcgd.cov <- function(y, x, var.type, lambda1 = NULL, lambda2 = NULL, U = NULL, maxit = 100, thresh = 1e-6,
                    alpha0 = NULL, theta0 = NULL, R0 = NULL, trace.it = T){
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y)
  p1 <- sum(var.type == "gaussian")
  p2 <- sum(var.type == "binomial")
  p3 <- sum(var.type == "poisson")
  groups <- as.factor(groups)
  if(is.null(alpha0)) {
    alpha0 <- rep(0, ncol(x))
    alpha.mat <-
      matrix(x%*% alpha0, nrow = n)
  } else{
    alpha.mat <-
      matrix(x%*% alpha0, nrow = n)
  }
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  param <- alpha.mat+theta0
  R0 <- svd(theta0, nu = 1, nv = 1)$d[1]
  if(is.null(lambda1)){
    sigma <- sum((y-mean(y, na.rm = T))^2, na.rm = T)/(n*p)
    beta <- sum(is.na(y))/sum(!is.na(y))*max(n,p)
    lambda1 <- 2*sigma*sqrt(beta*log(n+p))
  }
  low_bound <- 0
  ypois <- y0[, var.type == "poisson"]
  upper <- 12
  lower <- -12
  if(p2>0){
    low_bound <- low_bound - n*sum(var.type == "binomial")*(-2*max(abs(upper), abs(lower)) + log(1+exp(-2*max(abs(upper), abs(lower)))))
  }
  if(p3>0){
    low_bound <- low_bound - sum(ypois[!(ypois==0)]*(1-log(ypois[!(ypois==0)])), na.rm=T)
  }
  obj0 <- low_bound
  if(p1>0){
    obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                 na.rm = T)
  }
  if(p2>0){
    obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                         log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

  }
  if(p3>0){
    obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                         exp(param[, var.type == "poisson"]), na.rm = T)

  }
  U <- obj0/lambda1
  if(is.null(lambda2)) lambda2 <- 24*sum((y-mean(y))^2, na.rm = T)/(n*p)*log(n+p)
  alpha <- alpha0
  alpha.mat <-
    matrix(x %*% alpha, nrow = n)
  theta <- theta0
  R <- R0
  objective <- NULL
  error <- 1
  iter <- 0
  y0 <- y
  y0[is.na(y0)] <- 0
  low_bound <- 0
  ypois <- y0[, var.type == "poisson"]
  if(p2>0){
    low_bound <- low_bound - n*sum(var.type == "binomial")*(-2*max(abs(upper), abs(lower)) + log(1+exp(-2*max(abs(upper), abs(lower)))))
  }
  if(p3>0){
    low_bound <- low_bound - sum(ypois[!(ypois==0)]*(1-log(ypois[!(ypois==0)])))
  }

  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    alpha.tmp <- alpha
    theta.tmp <- theta
    R.tmp <- R
    grad <- grad_alpha(y, x, alpha, theta, var.type)
    step <- 1
    flag <- TRUE
    param <- alpha.mat + theta
    obj0 <- low_bound
    if(p1>0){
      obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                   na.rm = T)
    }
    if(p2>0){
      obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                           log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                           exp(param[, var.type == "poisson"]), na.rm = T)

    }
    obj0 <- obj0 + lambda2*sum(abs(alpha))+ lambda1*R
    while(flag){
      step <- 0.5*step
      mat <- abs(alpha.tmp - step *grad_alpha) - lambda2 * step
      alpha <- sign(alpha.tmp - step *grad_alpha) * pmax(as.matrix(mat), 0)
      alpha.mat <- matrix(x%*%alpha, nrow = n)
      param <- alpha.mat + theta
      obj <- low_bound
      if(p1>0){
        obj <- obj + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                   na.rm = T)
      }
      if(p2>0){
        obj <- obj + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                           log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

      }
      if(p3>0){
        obj <- obj + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                           exp(param[, var.type == "poisson"]), na.rm = T)

      }
      obj <- obj + lambda2*sum(abs(alpha))+ lambda1*R
      flag <- obj > obj0 +  thresh * abs(obj0)

    }
    alpha.mat <- matrix(x%*%alpha, nrow=n)
    grad_theta <- matrix(rep(0, n*p), nrow = n)
    if(p1>0){
      grad_theta[, var.type == "gaussian"] <- -omega[, var.type == "gaussian"]*(y0 - alpha.mat - theta)[, var.type == "gaussian"]
    }
    if(p2>0){
      grad_theta[, var.type == "binomial"] <- -omega[, var.type == "binomial"] *( y0 - exp(alpha.mat + theta)/(1+exp(alpha.mat + theta)))[, var.type == "binomial"]

    }
    if(p3>0){
      grad_theta[, var.type == "poisson"] <- - omega[, var.type == "poisson"]*(y0 - exp(alpha.mat + theta))[, var.type == "poisson"]

    }

    svd_theta <- rARPACK::svds(grad_theta, k=1, nu = 1, nv = 1)
    D_t <- - svd_theta$u%*%t(svd_theta$v)
    step <- 2
    flag <- TRUE
    param <- alpha.mat + theta
    obj0 <- low_bound
    if(p1>0){
      obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                   na.rm = T)
    }
    if(p2>0){
      obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                           log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                           exp(param[, var.type == "poisson"]), na.rm = T)

    }
    obj0 <- obj0 + lambda2*sum(abs(alpha))+ lambda1*R
    while((flag )){
      step <- 0.5*step
      if(lambda1 >= - sum(D_t*grad_theta)){
        R_hat <- 0
        theta_hat <- matrix(rep(0, n*p), nrow = n)
        if(norm(theta_hat - theta.tmp, type = "F")^2==0){
          beta <- 1
        } else {
          beta <- min(1, step)
        }
        theta <- theta.tmp + beta*(theta_hat - theta.tmp)
        R <- R.tmp + beta*(R_hat - R.tmp)
      } else{
        R_hat <- U
        theta_hat <- U*D_t
        if(norm(theta_hat - theta.tmp, type = "F")^2==0){
          beta <- 1
          theta <- theta.tmp + beta*(theta_hat - theta.tmp)
          R <- R.tmp + beta*(R_hat - R.tmp)
        } else {
          beta <- min(1, step)
          theta <- theta.tmp + beta*(theta_hat - theta.tmp)
          R <- R.tmp + beta*(R_hat - R.tmp)
        }
      }
      param <- alpha.mat + theta
      obj <- low_bound
      if(p1>0){
        obj <- obj + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                   na.rm = T)
      }
      if(p2>0){
        obj <- obj + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                           log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

      }
      if(p3>0){
        obj <- obj + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                           exp(param[, var.type == "poisson"]), na.rm = T)

      }
      obj <- obj + lambda2*sum(abs(alpha))+ lambda1*R
      flag <- (obj > obj0 + thresh * abs(obj0))
      if(step<=1e-10)  flagflag <- TRUE
    }
    #
    param <- alpha.mat + theta
    obj <- low_bound
    if(p1>0){
      obj <- obj + (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                                 na.rm = T)
    }
    if(p2>0){
      obj <- obj + sum(- (y[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                         log(1 + exp(param[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj <- obj + sum(- (y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                         exp(param[, var.type == "poisson"]), na.rm = T)

    }
    obj <- obj + lambda2*sum(abs(alpha))+ lambda1*R
    objective <- c(objective, obj)
    U <- obj/lambda1
    if(iter == 1){
      error <- 1
    } else {
      if(abs(objective[iter]) == 0) {
        error <- abs(objective[iter]-objective[iter - 1])
      } else {
        error <- abs(objective[iter]-objective[iter - 1]) /abs(objective[iter])
      }
    }
    if(trace.it && (iter%%10 == 0) ){
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  param <- alpha.mat + theta
  yy <- param
  yy[, var.type == "poisson"] <- exp(yy[, var.type == "poisson"])
  yy[, var.type == "binomial"] <- exp(yy[, var.type == "binomial"])/(1+exp(yy[, var.type == "binomial"]))
  y.imputed <- y
  y.imputed[is.na(y)] <- yy[is.na(y)]

  return(list(y.imputed = y.imputed,
              param=param,
              alpha = alpha,
              theta = theta))
}

