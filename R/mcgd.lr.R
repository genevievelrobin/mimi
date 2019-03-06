mcgd.lr <- function(y,
                    var.type,
                    lambda1 = NULL,
                    maxit = 100,
                    thresh = 1e-6,
                    theta0 = NULL,
                    trace.it = T){
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  p1 <- sum(var.type == "gaussian")
  p2 <- sum(var.type == "binomial")
  p3 <- sum(var.type == "poisson")
  omega <- !is.na(y)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  R0 <- svd(theta0, nu = 1, nv = 1)$d[1]
  if(is.null(lambda1)){
    sigma <- sum((y-mean(y, na.rm = T))^2, na.rm = T)/(n*p)
    beta <- sum(is.na(y))/sum(!is.na(y))*max(n,p)
    lambda1 <- 2*sigma*sqrt(beta*log(n+p))
  }
  theta <- theta0
  R <- R0
  objective <- NULL
  error <- 1
  iter <- 0
  list.theta <- list()
  list.R <- list()
  y0 <- y
  y0[is.na(y0)] <- 0
  low_bound <- 0
  ypois <- y0[, var.type == "poisson"]
  upper <- 12
  lower <- -12
  if(p2>0){
    low_bound <- low_bound - n*sum(var.type == "binomial")*(-2*max(abs(upper), abs(lower)) + log(1+exp(-2*max(abs(upper), abs(lower)))))
  }
  if(p3>0){
    low_bound <- low_bound - sum(ypois[!(ypois==0)]*(1-log(ypois[!(ypois==0)])))
  }
  obj0 <- low_bound
  if(p1>0){
    obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                                 na.rm = T)
  }
  if(p2>0){
    obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * theta[, var.type == "binomial"]) +
                         log(1 + exp(theta[, var.type == "binomial"])), na.rm = T)

  }
  if(p3>0){
    obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                         exp(theta[, var.type == "poisson"]), na.rm = T)

  }
  U <- obj0/lambda1
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    list.theta[[iter]] <- theta
    list.R[[iter]] <- R
    theta.tmp <- theta
    R.tmp <- R
    grad <- matrix(rep(0, n*p), nrow = n)
    if(p1>0){
      grad[, var.type == "gaussian"] <- -omega[, var.type == "gaussian"]*(y0 - theta)[, var.type == "gaussian"]
    }
    if(p2>0){
      grad[, var.type == "binomial"] <- -omega[, var.type == "binomial"] *( y0 - exp(theta)/(1+exp(theta)))[, var.type == "binomial"]

    }
    if(p3>0){
      grad[, var.type == "poisson"] <- - omega[, var.type == "poisson"]*(y0 - exp(theta))[, var.type == "poisson"]

    }
    step <- 2
    flag <- TRUE
    obj0 <- low_bound
    if(p1>0){
      obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                                   na.rm = T)
    }
    if(p2>0){
      obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * theta[, var.type == "binomial"]) +
                           log(1 + exp(theta[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                           exp(theta[, var.type == "poisson"]), na.rm = T)

    }
    obj0 <- obj0 + lambda1*R
    grad_theta <- matrix(rep(0, n*p), nrow = n)
    if(p1>0){
      grad_theta[, var.type == "gaussian"] <- -omega[, var.type == "gaussian"]*(y0 - theta)[, var.type == "gaussian"]
    }
    if(p2>0){
      grad_theta[, var.type == "binomial"] <- -omega[, var.type == "binomial"] *( y0 - exp(theta)/(1+exp(theta)))[, var.type == "binomial"]

    }
    if(p3>0){
      grad_theta[, var.type == "poisson"] <- - omega[, var.type == "poisson"]*(y0 - exp(theta))[, var.type == "poisson"]

    }

    svd_theta <- rARPACK::svds(grad_theta, k=1, nu = 1, nv = 1)
    D_t <- - svd_theta$u%*%t(svd_theta$v)
    step <- 2
    flag <- TRUE
    obj0 <- low_bound
    if(p1>0){
      obj0 <- obj0 + (1 / 2) * sum((y[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                                   na.rm = T)
    }
    if(p2>0){
      obj0 <- obj0 + sum(- (y[, var.type == "binomial"] * theta[, var.type == "binomial"]) +
                           log(1 + exp(theta[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj0 <- obj0 + sum(- (y[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                           exp(theta[, var.type == "poisson"]), na.rm = T)

    }
    obj0 <- obj0 + lambda1*R
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
      obj <- low_bound
      if(p1>0){
        obj <- obj + (1 / 2) * sum((y[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                                   na.rm = T)
      }
      if(p2>0){
        obj <- obj + sum(- (y[, var.type == "binomial"] * theta[, var.type == "binomial"]) +
                           log(1 + exp(theta[, var.type == "binomial"])), na.rm = T)

      }
      if(p3>0){
        obj <- obj + sum(- (y[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                           exp(theta[, var.type == "poisson"]), na.rm = T)

      }
      obj <- obj + lambda1*R
      flag <- (obj > obj0 + thresh * abs(obj0))
    }
    #
    obj <- low_bound
    if(p1>0){
      obj <- obj + (1 / 2) * sum((y[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                                 na.rm = T)
    }
    if(p2>0){
      obj <- obj + sum(- (y[, var.type == "binomial"] * theta[, var.type == "binomial"]) +
                         log(1 + exp(theta[, var.type == "binomial"])), na.rm = T)

    }
    if(p3>0){
      obj <- obj + sum(- (y[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                         exp(theta[, var.type == "poisson"]), na.rm = T)

    }
    obj <- obj + lambda1*R
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
  yy <- theta
  yy[, var.type == "poisson"] <- exp(yy[, var.type == "poisson"])
  yy[, var.type == "binomial"] <- exp(yy[, var.type == "binomial"])/(1+exp(yy[, var.type == "binomial"]))
  y.imputed <- y
  y.imputed[is.na(y)] <- yy[is.na(y)]

  return(list(y.imputed = y.imputed,
              theta = theta))
}

