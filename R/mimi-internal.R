#' irwls.lr
#'
#' @param y nxp observation matrix
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param nu positive number, backtracking line search parameter, default 0.01

#' @return A list vith the following elements
#' \item{y.imputed}{the imputed data set}
#' \item{theta}{estimated parameter matrix}
#' @import glmnet
irwls.lr <-
  function(y,
           var.type,
           lambda1,
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-5,
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
          th=0.1
        )
      theta <- res$theta
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
    estim[, var.type == "poisson"] <-
      matrix(rpois(n*sum(var.type == "poisson"),lambda=c(exp(estim[, var.type == "poisson"]))), nrow=n)
    estim[, var.type == "binary"] <- round(exp(estim[, var.type == "binary"]) / (1 + exp(estim[, var.type == "binary"])))
    y.imputed <- y0
    y.imputed[omega==0] <- estim[omega==0]
    return(list(
      theta = theta,
      y.estim = estim,
      y.imputed = y.imputed
    ))
  }

#' armijo.lr
#' Performs Armijo backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param w weights of the quadratic approximation
#' @param z matrix around which the quadratic approximation is done
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @param nu positive number, backtracking line search parameter, default 0.01
#' @param zeta positive number, backtracking line search parameter, default 0.1
#' @param th positive number, backtracking line search parameter, default 0.1
#' @import stats
#' @return A list with the following elements
#' \item{theta}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
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
           th=0.1) {
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
      sum(-(y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
            log(1 + exp(param.tmp[, var.type == "binary"])), na.rm = T)
    d.tmp <- svd(theta.tmp)$d
    flag <- TRUE
    step <- 1
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
        sum(-(y0[, var.type == "binary"] * param[, var.type == "binary"]) +
              log(1 + exp(param[, var.type == "binary"])), na.rm = T)
      d <- svd(theta)$d
      diff <- pois-pois.tmp
      diff <- diff + gaus-gaus.tmp
      diff <- diff + binom- binom.tmp
      diff <- diff + lambda1 * (sum(d) - sum(d.tmp))
      if (diff <= step * zeta * (-2 * sum(omega*w * z * direction, na.rm = T) +
                                th * sum((w*omega + (1-omega)*nu) * direction ^ 2) + lambda1 * (sum(d0) - sum(d.tmp)))) {
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

#' irwls.cov
#'
#' @param y nxp observation matrix
#' @param x (np)xN covariates matrix
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param alpha0 matrix of size (nb of groups)x(number of variables), initial value of the group effect, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param nu positive number, backtracking line search parameter, default 0.01

#' @return A list vith the following elements
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' @import glmnet
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
           nu=0.1) {
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n,p)-1)
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
    gaus <- (1 / 2) * sum(y0[, var.type == "gaussian"]^2, na.rm = T)
    pois <- n*p
    binom <-binom <--log(2)*n*p
    objective <- gaus+pois+binom
    y0 <- y
    error <- 1
    iter <- 0
    while (((error > thresh) && (iter < maxit))) {
      iter <- iter + 1
      alpha.tmp <- alpha
      theta.tmp <- theta
      yv <- quad_approx(y0, param, var.type)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega == 0] <- 0
      alpha <-
        glmnet(
          x,
          c((ytilde - omega*theta.tmp) * vtilde2 / (omega*vtilde2 + (1-omega)*nu)),
          family = "gaussian",
          lambda = lambda2,
          intercept = FALSE,
          weights = c(omega*vtilde2 + (1-omega)*nu)
        )

      alpha <- as.numeric(alpha$beta)
      alpha <-
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
          nu=nu,
          zeta=0.1,
          th=0.1
        )$alpha
      alpha.mat <-
        matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
      alpha.tmp <- alpha
      param <- alpha.mat + theta.tmp
      yv <- quad_approx(y0, param, var.type)
      ytilde <- yv$ytilde
      vtilde2 <- yv$vtilde2
      ytilde[omega == 0] <- 0
      w <- omega*vtilde2 + (1-omega)*nu
      lambda1w <- lambda1 / max(w)
      w <- w / max(w)
      svd_theta <-
        wlra(vtilde2 / (omega*vtilde2 + (1-omega)*nu) * (ytilde - omega * alpha.mat),
             w, lambda1w)
      u <- svd_theta$u
      d <- svd_theta$d
      v <- svd_theta$v
      if (is.null(dim(u))) {
        theta <- d * u %*% t(v)
      } else {
        theta <- u %*% diag(d) %*% t(v)
      }
      theta <-
        armijo.lr(y0=y,
                  theta=theta,
                  theta.tmp=theta.tmp,
                  alpha.mat=alpha.mat,
                  w=vtilde2,
                  z=ytilde - omega * param,
                  lambda1=lambda1,
                  var.type=var.type,
                  b=0.5,
                  nu=nu,
                  zeta=0.1,
                  th=0.1)$theta
      param <- alpha.mat + theta
      gaus <-
        (1 / 2) * sum((y[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                      na.rm = T)
      pois <-
        sum((-(y[, var.type == "poisson"] * param[, var.type == "poisson"]) +
               exp(param[, var.type == "poisson"])), na.rm = T)
      binom <-
        sum((-(y[, var.type == "binary"] * param[, var.type == "binary"]) +
               log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
      d <- svd(theta)$d
      objective <- c(objective, min(
        .Machine$double.xmax,
        (pois + gaus + binom + lambda1 * sum(d) + lambda2 * sum(abs(alpha)))
      ))
      if (iter == 1) {
        error <- 1
      } else{
        if (objective[iter] == .Machine$double.xmax) {
          error <- 1
        } else
          error <-
            abs(objective[iter] - objective[iter - 1]) / abs(objective[iter])

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
      matrix(rpois(n*sum(var.type == "poisson"),lambda=c(exp(estim[, var.type == "poisson"]))), nrow=n)
    estim[, var.type == "binary"] <- round(exp(estim[, var.type == "binary"]) / (1 + exp(estim[, var.type == "binary"])))
    y.imputed <- y0
    y.imputed[omega==0] <- estim[omega==0]
    return(
      list(
        y.imputed = y.imputed,
        param = alpha.mat + theta,
        alpha = alpha,
        theta = theta
      )
    )
  }


#' armijo.alpha
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param x (np)xN matrix of covariates
#' @param alpha  direction of update for vector of regression parameters of length N
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param alpha.tmp length N vector, current regression parameters
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param w weights of the quadratic approximation
#' @param z matrix around which the quadratic approximation is done
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda2 positive number, regularization parameter for l1 norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @param nu positive number, backtracking line search parameter, default 0.01
#' @param zeta positive number, backtracking line search parameter, default 0.1
#' @param th positive number, backtracking line search parameter, default 0.1
#' @import stats
#' @return A list with the following elements
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta a}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
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
           th=0.9) {
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
      sum((-(y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
             log(1 + exp(param.tmp[, var.type == "binary"]))), na.rm = T)
    flag <- TRUE
    a0 <- sum(abs(alpha.tmp + direction))
    step <- 1
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
        sum((-(y0[, var.type == "binary"] * param[, var.type == "binary"]) +
               log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
      diff <- gaus - gaus.tmp
      diff <- diff + pois - pois.tmp
      diff <- diff + binom - binom.tmp
      diff <- diff + lambda2 * (sum(abs(alpha)) - sum(abs(alpha.tmp)))
      if (diff <= step * zeta * (-2 * sum(omega*w * z * alpha.tmp.mat) + th *
                                sum((omega*w + (1-omega)*nu) * alpha.tmp.mat ^ 2) + lambda2 * (a0 - sum(abs(alpha.tmp))))) {
        flag <- FALSE
      } else
        step <- b * step
    }
    return(list(alpha = alpha))
  }


#' wght
#'
#' @param y real number: observation
#' @param param real number: current value of parameter
#' @param var.type type of variable y (gaussian, binary, poisson)
#' @return weight of the quadratic approximation
wght <- function(y, param, var.type) {
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
#' @return matrix of weights for the quadratic approximation (nxP)
quad_approx <- function(y, param, var.type) {
  n <- nrow(y)
  p <- length(var.type)
  yt <- rep(0, n)
  vt <- rep(0, n)
  for (j in 1:p) {
    w <-
      lapply(1:n, function(i)
        wght(
          y = y[i, j],
          param = param[i, j],
          var.type = var.type[j]
        ))
    ytilde <- sapply(1:n, function(i)
      w[[i]]$ytilde)
    vtilde2 <- sapply(1:n, function(i)
      w[[i]]$vtilde2)
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

