#' wls.cov
#'
#' @param y observation matrix
#' @param x covariates matrix
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param weights matrix of size (nb of ind.)x(number of variables) with entries in (0,1), optional weights matrix
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param mu0 real number, initial value of the offset, default 0
#' @param alpha0 a vector of length (nc covariates): initial regression parameter, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayed, default FALSE
#' @param offset boolean, if TRUE offset is computed, otherwise set to 0, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix
#'
#' @return A list vith the following elements
#' \item{y}{the original data matrix}
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a vector of length (nc covariates): the regression parameter}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' @import glmnet corpcor
wls.cov <-
  function(y,
           x,
           lambda1,
           lambda2,
           weights = NULL,
           thresh = 1e-5,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           trace.it = F,
           offset = F,
           max.rank = 5) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    q <- ncol(x)
    y0 <- y
    omega <- !is.na(y)
    if (is.null(mu0))
      mu0 <- 0
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    mu <- mu0
    alpha <- alpha0
    theta <- theta0
    objective <- NULL
    y <- y0
    y[!omega] <- 0
    error <- 1
    iter <- 0
    if (is.null(weights))
      weights <- matrix(rep(1, n * p), nrow = n)
    while ((error > thresh) && (iter < maxit)) {
      iter <- iter + 1
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      y.tmp <- mu + alpha.mat + theta
      y <- y0
      y[!omega] <- y.tmp[!omega]
      mu.tmp <- mu
      alpha.tmp <- alpha
      theta.tmp <- theta
      alpha <-
        glmnet(
          x,
          c(y - mu - theta.tmp),
          family = "gaussian",
          lambda = lambda2,
          intercept = FALSE,
          thresh = thresh,
          weights = c(weights)
        )
      mu <- alpha$a0
      alpha <- as.numeric(alpha$beta)
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      svd_theta <-
        wlra(
          y - mu - alpha.mat,
          w = weights,
          lambda = lambda1,
          x0 = NULL,
          thresh = thresh,
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
      objective <-
        c(objective,
          sum((1 / 2) * weights * (y0 - mu - alpha.mat - theta) ^ 2, na.rm = T) + lambda1 * sum(d) +
            lambda2 * sum(abs(alpha)))
      if (iter == 1) {
        error <- 1
      } else{
        error <-
          abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])
      }
      if (trace.it && (iter %% 10 == 0)) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    y <- mu + alpha.mat + theta
    y.imputed <- y0
    y.imputed[is.na(y0)] <- y[is.na(y0)]
    return(
      list(
        y = y0,
        y.imputed = y.imputed,
        param = mu + alpha.mat + theta,
        mu = mu,
        alpha = alpha,
        theta = theta,
        objective = objective,
        iter = iter
      )
    )

  }

#' irwls.cov
#'
#' @param y nxp observation matrix
#' @param x (np)xN covariates matrix
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param mu0 real number, initial value of the offset, default 0
#' @param alpha0 matrix of size (nb of groups)x(number of variables), initial value of the group effect, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param offset boolean, if TRUE offset is computed, otherwise set to 0, default FALSE
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)

#' @return A list vith the following elements
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' @import FactoMineR corpcor
irwls.cov <-
  function(y,
           x,
           var.type,
           lambda1,
           lambda2,
           nlevel = NULL,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           offset = F,
           scale = F,
           max.rank = 5,
           vt2 = NULL) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    wmax = 2 * max(y, na.rm = T)
    q <- ncol(x)
    x <- as.matrix(x)
    x <- matrix(as.numeric(x), nrow = nrow(x))
    omega <- !is.na(y)
    if (is.null(mu0))
      mu0 <- 0
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    mu <- mu0
    alpha <- alpha0
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    theta <- theta0
    param <- mu + alpha.mat + theta
    objective <- NULL
    y0 <- y
    error <- 1
    iter <- 0
    if (scale == T) {
      moy <- colMeans(y, na.rm = T)
      moy[vt2 == "binary"] <-
        log(moy[vt2 == "binary"] / (1 - moy[vt2 == "binary"]))
      moy[vt2 == "poisson"] <- log(moy[vt2 == "poisson"])
      moy[vt2 == "categorical"] <-
        log(moy[vt2 == "categorical"] / (1 - moy[vt2 == "categorical"]))
      scaling <- rep(0, p)
      p1 <- sum(var.type == "gaussian")
      p2 <- sum(var.type == "binary")
      p3 <- sum(var.type == "poisson")
      p4 <- sum(var.type == "categorical")
      if (p1 > 0) {
        scaling[vt2 == "gaussian"] <-
          0.5 * colSums((y[, vt2 == "gaussian"] - t(matrix(
            rep(moy[vt2 == "gaussian"],
                n), nrow = p1
          ))) ^ 2, na.rm = T) / (n - 1)
      }
      if (p2 > 0) {
        scaling[vt2 == "binary"] <-
          colSums(-(y[, vt2 == "binary"] * t(matrix(
            rep(moy[vt2 == "binary"],
                n), nrow = p2
          ))) +
            log(1 + exp(t(
              matrix(rep(moy[vt2 == "binary"],
                         n), nrow = p2)
            ))), na.rm = T) / (n - 1)

      }
      if (p3 > 0) {
        scaling[vt2 == "poisson"] <-
          colSums(-(y[, vt2 == "poisson"] * t(matrix(
            rep(moy[vt2 == "poisson"],
                n), nrow = p3
          ))) +
            exp(t(matrix(
              rep(moy[vt2 == "poisson"], n), nrow = p3
            ))) + log_factorial(y[, vt2 == "poisson"]), na.rm = T) / (n - 1)

      }
      if (p4 > 0) {
        scaling[vt2 == "categorical"] <-
          colSums(-(y[, vt2 == "categorical"] * t(matrix(
            rep(moy[vt2 == "categorical"],
                n), nrow = sum(vt2 == "categorical")
          ))) +
            log(1 + exp(t(
              matrix(rep(moy[vt2 == "categorical"],
                         n), nrow = sum(vt2 == "categorical"))
            ))), na.rm = T) / (n - 1)

      }
    } else
      scaling <- rep(1, ncol(y))
    sc <- matrix(rep(scaling, n), nrow = n, byrow = T)
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      mu.tmp <- mu
      alpha.tmp <- alpha
      theta.tmp <- theta
      yv <- quad_approx(y0, param, var.type, nlevel, wmax)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega=0] <- 0
      lambda1w <- lambda1 / (max(vtilde2))
      lambda2w <- lambda2 / (max(vtilde2))
      vtilde2 <- vtilde2 / max(vtilde2)
      if (scale == T) {
        vtilde2 <- sweep(vtilde2, 2, scaling, "/")
      }
      res_approx <-
        wls.cov(
          ytilde,
          x,
          lambda1w,
          lambda2w,
          weights = vtilde2,
          thresh = thresh,
          mu0 = mu.tmp,
          alpha0 = alpha.tmp,
          theta0 = theta.tmp,
          trace.it = F,
          offset = offset
        )
      mu <- res_approx$mu
      alpha <- res_approx$alpha
      theta <- res_approx$theta
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      param <- mu + alpha.mat + theta
      res <-
        bls.cov(
          y0,
          x,
          mu,
          alpha,
          theta,
          mu.tmp,
          alpha.tmp,
          theta.tmp,
          b = 0.5,
          lambda1,
          lambda2,
          var.type,
          thresh,
          sc = sc,
          nlevel = nlevel,
          vt2 = vt2
        )
      mu <- res$mu
      alpha <- res$alpha
      theta <- res$theta
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      param <- mu + alpha.mat + theta
      gaus <-
        (1 / 2) * sum(sc[, var.type == "gaussian"] * (y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum(sc[, var.type == "poisson"] * (-(y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                                             exp(param[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum(sc[, var.type == "binary"] * (-(y0[, var.type == "binary"] * param[, var.type == "binary"]) +
                                            log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
      truc <- rep(0, n)
      if (sum(var.type == "categorical") > 0) {
        for (j in 1:sum(var.type == "categorical")) {
          tt <-
            rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                             "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
          truc <-
            cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
        }
        truc <- truc[, 2:ncol(truc)]
        cat <-
          sum(sc[, vt2 == "categorical"] * (-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                log(truc))), na.rm = T)
      } else
        cat <- 0

      d <- svd(theta)$d
      if (length(lambda2) == 1) {
        objective <- c(objective, min(
          .Machine$double.xmax,
          (pois + gaus + binom + cat + lambda1 * sum(d) + lambda2 * sum(abs(alpha)))
        ))
      } else {
        objective <- c(objective, min(
          .Machine$double.xmax,
          (pois + gaus + binom + cat + lambda1 * sum(d) + sum(lambda2 * t(abs(
            alpha
          ))))
        ))
      }

      if (iter == 1) {
        error <- 1
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <-
            abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])

      }
      if (all(var.type == "gaussian"))
        error <- 0
      if (trace.it) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    y <- mu + alpha.mat + theta
    y[, var.type == "poisson"] <- exp(y[, var.type == "poisson"])
    y[, var.type == "binary"] <-
      exp(y[, var.type == "binary"]) / (1 + exp(y[, var.type == "binary"]))
    y.imputed <- y0
    y.imputed[is.na(y0)] <- y[is.na(y0)]
    return(
      list(
        y.imputed = y.imputed,
        param = mu + alpha.mat + theta,
        mu = mu,
        alpha = alpha,
        theta = theta
      )
    )
  }


#' irwls.lr
#'
#' @param y nxp observation matrix
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param offset boolean, if TRUE offset is computed, otherwise set to 0, default FALSE
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)

#' @return A list vith the following elements
#' \item{y.imputed}{the imputed data set}
#' \item{theta}{estimated parameter matrix}
#' @import FactoMineR corpcor
irwls.lr <-
  function(y,
           var.type,
           lambda1,
           nlevel = NULL,
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           offset = F,
           scale = F,
           max.rank = 5,
           vt2 = NULL) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    theta <- theta0
    param <- theta
    objective <- NULL
    y0 <- y
    error <- 1
    iter <- 0
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      theta.tmp <- theta
      param.tmp <- param
      yv <- quad_approx(y0, param, var.type, nlevel)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      if (scale == T) {
        scaling <- apply(y, 2, sd, na.rm = T)
        ytilde <- sweep(ytilde, 2, scaling, "/")
      }
      ytilde[omega==0] <- 0
      z <- ytilde*vtilde2/(vtilde2+0.05)
      w <- vtilde2+0.05
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
      if (sum(var.type == "categorical"))
        theta[, vt2 == "categorical"] <-
        sweep(theta[, vt2 == "categorical"], 1, rowMeans(theta[, vt2 == "categorical"]))
      if (scale) {
        theta <- sweep(theta, 2, scaling, "*")
      }
      res <-
        armijo.lr(y0, theta, theta.tmp, b = 0.5, lambda1, vt2, thresh, nlevel, vt2,
                  w=vtilde2, z=ytilde)
      theta <- res$theta
      param <- theta
      objective <-
        c(objective, min(.Machine$double.xmax, res$objective))
      if (iter == 1) {
        error <- 1
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <-
            abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])

      }
      if (all(var.type == "gaussian"))
        error <- 0
      if (trace.it) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    y <- theta
    y[, vt2 == "poisson"] <- exp(y[, vt2 == "poisson"])
    y[, vt2 == "binary"] <-
      exp(y[, vt2 == "binary"]) / (1 + exp(y[, vt2 == "binary"]))
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        count <- sum(nlevel[1:(which(var.type == "categorical")[j] - 1)])
        y[, (count + 1):(count + nlevel[which(vt2 == "categorical")[j]])] <-
          t(sapply(1:n, function(i)
            exp(y[i, (count + 1):(count + nlevel[which(vt2 == "categorical")[j]])]) /
              sum(exp(y[i, (count + 1):(count + nlevel[which(vt2 == "categorical")[j]])]))))
      }
    }
    y.imputed <- y0
    y.imputed[is.na(y.imputed)] <- y[is.na(y.imputed)]
    return(list(y.imputed = y.imputed, theta = theta))
  }



#' wls.multi
#'
#' @param y observation matrix
#' @param groups factor indicating group memberships, if vector will be treated as a factor
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param weights nxp matrix of weights in (0,1)
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param mu0 real number, initial value of the offset, default 0
#' @param alpha0 matrix of size (nb of groups)x(number of variables), initial value of the group effect, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param offset boolean, whether an offset should be fitted
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @import corpcor
#'
#' @return A list with the following elements
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta a}{(nb individuals) x (nb variables) matrix containing the individual effects}
wls.multi <- function(y,
                      groups,
                      lambda1,
                      lambda2,
                      weights = NULL,
                      thresh = 1e-5,
                      maxit = 1e3,
                      mu0 = NULL,
                      alpha0 = NULL,
                      theta0 = NULL,
                      trace.it = F,
                      offset = F,
                      max.rank = 5) {
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  y <- matrix(as.numeric(as.character(y)), nrow = n)
  y0 <- y
  wmax = 2 * max(y, na.rm = T)
  groups <- as.factor(groups)
  N <- nlevels(groups)
  omega <- !is.na(y)
  ncenters <- aggregate(rep(1, n), list(groups), sum)[, 2]
  if (is.null(weights))
    weights <- matrix(rep(1, n * p), nrow = n)
  if (is.null(mu0))
    mu0 <- 0
  if (is.null(alpha0)) {
    alpha0 <- rep(0, N * p)
    alpha.rep <-
      matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
    alpha0 <- matrix(as.matrix(alpha0), nrow = N)
  } else{
    alpha.rep <-
      matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
    alpha <- matrix(as.matrix(alpha0), nrow = N)
  }
  if (is.null(theta0))
    theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  alpha.rep <-
    matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
  theta <- theta0
  objective <- NULL
  w.groups <-
    aggregate(matrix(weights * rep(1, n * p), nrow = n), list(groups), sum)[, 2:(p + 1)]
  w.groups[w.groups == 0] <- 1
  y <- y0
  error <- 1
  iter <- 0
  while (((error > thresh) && (iter < maxit))) {
    iter <- iter + 1
    y.tmp <- mu + alpha.rep + theta
    y <- y0
    y[!omega] <- y.tmp[!omega]
    mu.tmp <- mu
    alpha.tmp <- alpha
    theta.tmp <- theta
    if (offset)
      mu <-
      colMeans(weights * (y - alpha.rep - theta.tmp), na.rm = T)
    else
      mu <- 0
    y.center <-
      aggregate(weights * (y - mu - theta.tmp), list(groups), sum, na.rm = T)[, 2:(p +
                                                                                     1)] / w.groups
    mat <-
      abs(as.matrix(y.center)) - lambda2 / (2 * as.matrix(w.groups))
    alpha <- sign(as.matrix(y.center)) * pmax(mat, 0)
    alpha.rep <-
      matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
    rep.y.center <-
      matrix(rep(as.matrix(y.center), rep(ncenters, p)), nrow = n)
    svd_theta <-
      wlra(
        y - mu - alpha.rep,
        w = weights,
        lambda = lambda1,
        thresh = thresh,
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
    objective <-
      c(objective,
        sum((1 / 2) * weights * (y0 - t(
          matrix(rep(mu, n), nrow = p)
        ) - alpha.rep - theta) ^ 2, na.rm = T) + lambda1 * sum(d) +
          lambda2 * sum(abs(alpha)))
    if (iter == 1) {
      error <- 1
    } else{
      error <-
        abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])
    }
    if (trace.it && (iter %% 10 == 0)) {
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  y <- mu + alpha.rep + theta
  y.imputed <- y0
  y.imputed[is.na(y0)] <- y[is.na(y0)]
  return(
    list(
      y.imputed = y.imputed,
      param = t(matrix(rep(mu, n), nrow = p)) + alpha.rep + theta,
      mu = mu,
      alpha = alpha,
      theta = theta
    )
  )
}


#' irwls.multi
#'
#' @param y observation matrix
#' @param groups factor indicating group memberships, if vector will be treated as a factor
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param mu0 real number, initial value of the offset, default 0
#' @param alpha0 matrix of size (nb of groups)x(number of variables), initial value of the group effect, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param offset boolean, whether an offset should be fitted
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)
#'
#' @return A list vith the following elements
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' @import FactoMineR corpcor
irwls.multi <-
  function(y,
           groups,
           var.type,
           lambda1,
           lambda2,
           nlevel = NULL,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           offset = F,
           scale = F,
           max.rank = 5,
           vt2 = NULL) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    wmax = 2 * max(y, na.rm = T)
    groups <- as.factor(groups)
    N <- nlevels(groups)
    ncenters <- aggregate(rep(1, n), list(groups), sum)[, 2]
    if (is.null(mu0))
      mu0 <- 0
    if (is.null(alpha0)) {
      alpha0 <- rep(0, N * p)
      alpha.rep <-
        matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
      alpha0 <- matrix(as.matrix(alpha0), nrow = N)
    } else{
      alpha.rep <-
        matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
      alpha <- matrix(as.matrix(alpha0), nrow = N)
    }
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    mu <- mu0
    alpha <- alpha0
    alpha.rep <-
      matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
    theta <- theta0
    param <- mu + alpha.rep + theta
    objective <- NULL
    y0 <- y
    error <- 1
    iter <- 0
    omega <- !is.na(y)
    ncenters <- aggregate(rep(1, n), list(groups), sum)[, 2]

    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      mu.tmp <- mu
      alpha.tmp <- alpha
      theta.tmp <- theta
      yv <- quad_approx(y0, param, var.type, nlevel, wmax)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      if (scale == T) {
        scaling <- apply(y, 2, sd, na.rm = T)
        ytilde <- sweep(ytilde, 2, scaling, "/")
      }
      lambda1w <- lambda1 / (max(vtilde2))
      lambda2w <- lambda2 / (max(vtilde2))
      vtilde2 <- vtilde2 / max(vtilde2)
      res_approx <-
        wls.multi(
          ytilde,
          groups,
          lambda1w,
          lambda2w,
          weights = vtilde2,
          thresh = thresh,
          mu0 = mu.tmp,
          alpha0 = alpha.tmp,
          theta0 = theta.tmp,
          trace.it = F,
          offset = offset,
          max.rank = max.rank
        )
      mu <- res_approx$mu
      alpha <- res_approx$alpha
      theta <- res_approx$theta
      alpha.rep <-
        matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
      if (scale) {
        theta <- sweep(theta, 2, scaling, "*")
      }
      param <- mu + alpha.rep + theta

      res <-
        bls.multi(
          y0,
          groups,
          mu,
          alpha,
          theta,
          mu.tmp,
          alpha.tmp,
          theta.tmp,
          b = 0.5,
          lambda1,
          lambda2,
          var.type,
          thresh,
          nlevel,
          vt2
        )
      mu <- res$mu
      alpha <- res$alpha
      theta <- res$theta
      alpha.rep <-
        matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
      param <- mu + alpha.rep + theta
      objective <-
        c(objective, min(.Machine$double.xmax, res$objective))
      if (iter == 1) {
        error <- 1
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <-
            abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])

      }
      if (all(var.type == "gaussian"))
        error <- 0
      if (trace.it) {
        print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
      }
    }
    if (error < thresh)
      cvg = T
    else
      cvg = F
    y <- mu + alpha.rep + theta
    y[, var.type == "poisson"] <- exp(y[, var.type == "poisson"])
    y[, var.type == "binary"] <-
      exp(y[, var.type == "binary"]) / (1 + exp(y[, var.type == "binary"]))
    y.imputed <- y0
    y.imputed[is.na(y0)] <- y[is.na(y0)]
    return(
      list(
        y.imputed = y.imputed,
        param = mu + alpha.rep + theta,
        mu = mu,
        alpha = alpha,
        theta = theta
      )
    )
  }



#' wght
#'
#' @param y real number: observation
#' @param param real number: current value of parameter
#' @param var.type type of variable y (gaussian, binary, poisson)
#' @param nlevel vector indicating, for every column of y, the number of parameters of the distribution (1 for numerics, K for categorical with K categories)
#' @return weight of the quadratic approximation
wght <- function(y, param, var.type, nlevel = NULL) {
  if (var.type == "gaussian") {
    ytilde <- y
    vtilde2 <- 1
  } else if (var.type == "binary") {
    vtilde2 <- 1 / 4
    ytilde <-
      y / vtilde2 - (exp(param) / (1 + exp(param))) / vtilde2 + param
  } else if (var.type == "poisson") {
    vtilde2 <- exp(param)
    ytilde <- (y - exp(param)) / vtilde2  + param
  } else if (var.type == "categorical") {
    vtilde2 <- rep(0.25, nlevel)
    ytilde <-
      y / vtilde2 - (exp(param) / sum(exp(param))) / vtilde2 + param
  }
  else {
    print(var.type)
    stop(
      "Incorrect type of variable. Should be 'gaussian', 'binary', 'poisson' or 'categorical'."
    )
  }
  return(list(ytilde = ytilde, vtilde2 = vtilde2))
}

#' quad_approx
#'
#' @param y observation matrix (nxP) - extended with dummies for cat variables
#' @param param matrix of parameters (nxP) - extended for cat variables
#' @param var.type type of the variables in y (length p<=P) - NOT extended for cat variables
#' @param nlevel vector indicating, for every column of y, the number of parameters of the distribution (1 for numerics, K for categorical with K categories) - NOT extended for cat variables
#' @param wmax maximum weight of quadratic approximation
#' @return matrix of weights for the quadratic approximation (nxP)
quad_approx <- function(y, param, var.type, nlevel, wmax) {
  n <- nrow(y)
  p <- length(var.type)
  yt <- rep(0, n)
  vt <- rep(0, n)
  count <- 1
  for (j in 1:p) {
    w <-
      lapply(1:n, function(i)
        wght(
          y = y[i, count:(count + nlevel[j] - 1)],
          param = param[i, count:(count + nlevel[j] - 1)],
          var.type = var.type[j],
          nlevel = nlevel[j]
        ))
    count <- count + nlevel[j]
    if (nlevel[j] == 1) {
      ytilde <- sapply(1:n, function(i)
        w[[i]]$ytilde)
      vtilde2 <- sapply(1:n, function(i)
        w[[i]]$vtilde2)
    } else{
      ytilde <- t(sapply(1:n, function(i)
        w[[i]]$ytilde))
      vtilde2 <- t(sapply(1:n, function(i)
        w[[i]]$vtilde2))
    }
    yt <- cbind(yt, as.matrix(ytilde))
    vt <- cbind(vt, as.matrix(vtilde2))
  }
  return(list(ytilde = yt[, 2:ncol(yt)], vtilde2 = vt[, 2:ncol(vt)]))
}


#' log_factorial
#'
#' @param x matrix of integers
#' @return the log of the factorial of every entry in x
log_factorial <- function(x) {
  m <- apply(x, c(1, 2), function(t)
    if (is.na(t)) {
      NA
    } else if (t == 0) {
      1
    } else
      sum(log(1:t)))
  return(m)
}

#' wlra
#'
#' @param x nxp matrix to be approximates
#' @param w nxp matrix of weights (optional)
#' @param lambda regularization parameter
#' @param x0 nxp matrix of initial values
#' @param thresh convergence criterion
#' @param maxit maximum number of iterations
#' @param rank.max maximum number of sv to compute
#' @return a list containing the following elements
#' \item{d}{the vector of singular values}
#' \item{u}{the left singular vectors}
#' \item{v}{the right singular vectors}
#' \item{convergence}{boolean indicating if algorithm converged before maxit iterations}
#' \item{iter}{number of iterations performed}
#' @import corpcor
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
      svd.old = svd(xhat)
      svd.old$xhat.old <- xhat
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
    if (iter == maxit)
      warning(paste("Convergence not achieved by", maxit, "iterations"))
    return(list(
      d = svd.xfill$d,
      u = svd.xfill$u,
      v = svd.xfill$v,
      cvg = cvg,
      iter = iter
    ))
  }

#' bls.cov
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param x (np)xN matrix of covariates
#' @param mu real number, direction of update for offset
#' @param alpha  direction of update for vector of regression parameters of length N
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param mu.tmp real number, current offset
#' @param alpha.tmp length N vector, current regression parameters
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param lambda2 positive number, regularization parameter for l1 norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)
#' @param sc scaling matrix
#' @import stats corpcor
#' @return A list with the following elements
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta a}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
bls.cov <-
  function(y0,
           x,
           mu,
           alpha,
           theta,
           mu.tmp,
           alpha.tmp,
           theta.tmp,
           b = 0.5,
           lambda1,
           lambda2,
           var.type,
           thresh = 1e-5,
           sc,
           nlevel,
           vt2) {
    d <- dim(y0)
    n <- d[1]
    p <- d[2]
    q <- ncol(x)
    omega <- !is.na(y0)
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    alpha.tmp.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha.tmp, nrow = n)
    param <- mu + alpha.mat + theta
    param.tmp <- mu.tmp + alpha.tmp.mat + theta.tmp
    gaus.tmp <-
      (1 / 2) * sum(sc[, var.type == "gaussian"] * (y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois.tmp <-
      sum((-(y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
             exp(param.tmp[, var.type == "poisson"])), na.rm = T)
    binom.tmp <-
      sum((-(y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
             log(1 + exp(param.tmp[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <-
          rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                           "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
        truc <-
          cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat.tmp <-
        sum((-(y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                 log(truc))), na.rm = T)
    } else
      cat.tmp <- 0
    d.tmp <- svd(theta.tmp)$d
    gaus <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois <-
      sum((-(y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
             exp(param[, var.type == "poisson"])), na.rm = T)
    binom <-
      sum((-(y0[, var.type == "binary"] * param[, var.type == "binary"]) +
             log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <-
          rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                           "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
        truc <-
          cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <-
        sum((-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                 log(truc))), na.rm = T)
    } else
      cat <- 0
    d <- svd(theta)$d
    t <- 1
    mu2 <- (1 - t) * mu.tmp + t * mu
    alpha2 <- (1 - t) * alpha.tmp + t * alpha
    theta2 <- (1 - t) * theta.tmp + t * theta
    param2 <- (1 - t) * param.tmp + t * param
    diff <-
      pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha))))
    number <-
      gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp) + sum(lambda2 * (t(abs(alpha.tmp))))
    while (diff < -abs(number) * thresh) {
      t <- b * t
      mu2 <- (1 - t) * mu.tmp + t * mu
      alpha2 <- (1 - t) * alpha.tmp + t * alpha
      theta2 <- (1 - t) * theta.tmp + t * theta
      param2 <- (1 - t) * param.tmp + t * param
      gaus <-
        (1 / 2) * sum((y0[, var.type == "gaussian"] - param2[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum((-(y0[, var.type == "poisson"] * param2[, var.type == "poisson"]) +
               exp(param2[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum((-(y0[, var.type == "binary"] * param2[, var.type == "binary"]) +
               log(1 + exp(param2[, var.type == "binary"]))), na.rm = T)
      truc <- rep(0, n)
      if (sum(var.type == "categorical") > 0) {
        for (j in 1:sum(var.type == "categorical")) {
          tt <-
            rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                             "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
          truc <-
            cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
        }
        truc <- truc[, 2:ncol(truc)]
        cat <-
          sum((-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                   log(truc))), na.rm = T)
      } else
        cat <- 0
      d <- svd(theta2)$d
      diff <-
        pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha2))))
    }
    obj <-
      pois + gaus + binom + lambda1 * d + sum(lambda2 * t(abs(alpha2)))
    return(list(
      mu = mu2,
      alpha = alpha2,
      theta = theta2,
      objective = obj,
      t = t
    ))

  }

#' bls.multi
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp matrix of observations
#' @param groups length n indicator factor of group memberships
#' @param mu real number new offset
#' @param alpha vector of length p*N(nb of groups) new regression parameter
#' @param theta nxp matrix new interaction matrix
#' @param mu.tmp real number current offset
#' @param alpha.tmp vector of length p*N(nb of groups) current regression parameter
#' @param theta.tmp nxp matrix current interaction matrix
#' @param b number in (0,1) factor by which step size is reduced
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param var.type vector of length p indicating variable types (gaussian, binary, poisson)
#' @param thresh positive number congervence criterion
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)
#'
#' @import stats corpcor
bls.multi <-
  function(y0,
           groups,
           mu,
           alpha,
           theta,
           mu.tmp,
           alpha.tmp,
           theta.tmp,
           lambda1,
           lambda2,
           var.type,
           thresh = 1e-5,
           b = 0.5,
           nlevel,
           vt2) {
    d <- dim(y0)
    n <- d[1]
    p <- d[2]
    groups <- as.factor(groups)
    N <- nlevels(groups)
    omega <- !is.na(y0)
    ncenters <- aggregate(rep(1, n), by = list(groups), sum)[, 2]
    alpha.rep <-
      matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
    alpha.tmp.rep <-
      matrix(rep(as.matrix(alpha.tmp), rep(ncenters, p)), nrow = n)
    param <- mu + alpha.rep + theta
    param.tmp <- mu.tmp + alpha.tmp.rep + theta.tmp
    gaus.tmp <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois.tmp <-
      sum(-y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"] +
            exp(param.tmp[, var.type == "poisson"]), na.rm = T)
    binom.tmp <-
      sum(-y0[, var.type == "binary"] * param.tmp[, var.type == "binary"] +
            log(1 + exp(param.tmp[, var.type == "binary"])), na.rm = T)
    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <-
          rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                           "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
        truc <-
          cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat.tmp <-
        sum(-(y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                log(truc)), na.rm = T)
    } else
      cat.tmp <- 0
    d.tmp <- svd(theta.tmp)$d
    gaus <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois <-
      sum(-y0[, var.type == "poisson"] * param[, var.type == "poisson"] +
            exp(param[, var.type == "poisson"]), na.rm = T)
    binom <-
      sum(-y0[, var.type == "binary"] * param[, var.type == "binary"] +
            log(1 + exp(param[, var.type == "binary"])), na.rm = T)

    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <-
          rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                           "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
        truc <-
          cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <-
        sum((-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                 log(truc))), na.rm = T)
    } else
      cat <- 0
    d <- svd(theta)$d
    t <- 1
    mu2 <- (1 - t) * mu.tmp + t * mu
    alpha2 <- (1 - t) * alpha.tmp + t * alpha
    theta2 <- (1 - t) * theta.tmp + t * theta
    param2 <- (1 - t) * param.tmp + t * param
    diff <-
      pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha))))
    number <-
      gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp) + sum(lambda2 * (t(abs(alpha.tmp))))
    while (diff < -abs(number) * thresh / 2) {
      t <- b * t
      mu2 <- (1 - t) * mu.tmp + t * mu
      alpha2 <- (1 - t) * alpha.tmp + t * alpha
      theta2 <- (1 - t) * theta.tmp + t * theta
      param2 <- (1 - t) * param.tmp + t * param
      gaus <-
        (1 / 2) * sum((y0[, var.type == "gaussian"] - param2[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum((-(y0[, var.type == "poisson"] * param2[, var.type == "poisson"]) +
               exp(param2[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum((-(y0[, var.type == "binary"] * param2[, var.type == "binary"]) +
               log(1 + exp(param2[, var.type == "binary"]))), na.rm = T)
      truc <- rep(0, n)
      if (sum(var.type == "categorical") > 0) {
        for (j in 1:sum(var.type == "categorical")) {
          tt <-
            rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                             "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
          truc <-
            cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
        }
        truc <- truc[, 2:ncol(truc)]
        cat <-
          sum((-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                   log(truc))), na.rm = T)
      } else
        cat <- 0
      d <- svd(theta2)$d

      diff <-
        pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha2))))
    }
    obj <-
      pois + gaus + binom + lambda1 * d + sum(lambda2 * t(abs(alpha2)))
    return(list(
      mu = mu2,
      alpha = alpha2,
      theta = theta2,
      objective = obj,
      t = t
    ))
  }


#' armijo.lr
#' Performs Armijo backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)
#' @import stats corpcor
#' @return A list with the following elements
#' \item{theta}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
armijo.lr <-
  function(y0,
           theta,
           theta.tmp,
           w,
           z,
           b = 0.5,
           lambda1,
           var.type,
           thresh = 1e-5,
           nlevel,
           vt2) {
    d <- dim(y0)
    n <- d[1]
    p <- d[2]
    omega <- !is.na(y0)
    param <- theta
    param.tmp <- theta.tmp
    gaus.tmp <-
      (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"]) ^
                      2,
                    na.rm = T)
    pois.tmp <-
      sum(-(y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
            exp(param.tmp[, var.type == "poisson"]), na.rm = T)
    binom.tmp <-
      sum(-(y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
            log(1 + exp(param.tmp[, var.type == "binary"])), na.rm = T)
    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <-
          rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                           "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
        truc <-
          cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat.tmp <-
        sum(-(y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                log(truc)), na.rm = T)
    } else
      cat.tmp <- 0

    d.tmp <- svd(theta.tmp)$d
    flag <- TRUE
    step <- 1
    print(step)
    while(flag){
      print(step)
      param <- param.tmp + step*param
      gaus <-
        (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"]) ^
                        2,
                      na.rm = T)
      pois <-
        sum(-(y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
              exp(param[, var.type == "poisson"]), na.rm = T)
      binom <-
        sum(-(y0[, var.type == "binary"] * param[, var.type == "binary"]) +
              log(1 + exp(param[, var.type == "binary"])), na.rm = T)

      truc <- rep(0, n)
      if (sum(var.type == "categorical") > 0) {
        for (j in 1:sum(var.type == "categorical")) {
          tt <-
            rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type ==
                                                                             "categorical")[j] + nlevel[which(var.type == "categorical")[j]] - 1)]))
          truc <-
            cbind(truc, matrix(rep(tt, nlevel[which(var.type == "categorical")[j]]), nrow = n))
        }
        truc <- truc[, 2:ncol(truc)]
        cat <-
          sum(-(y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                  log(truc)), na.rm = T)
      } else
        cat <- 0
      d <- svd(param)$d
      diff <- pois+gaus+binom+cat-pois.tmp-gaus.tmp-binom.tmp-cat.tmp+lambda1*(sum(d)-sum(d.tmp))
      if(diff<=-2*sum(w*z*param.tmp, na.rm=T)+0.1*sum((w+0.05)*param.tmp^2)+lambda1*(sum(d)-sum(d.tmp))){
        flag <- FALSE
      } else{
        step <- step/2
      }
    }
    obj <- pois + gaus + binom + lambda1 * d
    return(list(
      theta = param,
      objective = obj,
      step = step
    ))

  }
