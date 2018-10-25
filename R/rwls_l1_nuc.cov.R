#' wls_l1_nuc.cov
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
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
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
#' @export
#' @import glmnet
#'
#' @examples
#' y <- matrix(rnorm(6 * 10), nrow = 6)
#' y[sample(1:50, size = 10)] <- NA
#' x <- matrix(rnorm(60*2), nrow = 60)
#' res <- wls_l1_nuc.cov(y, x, lambda1 = 0.1, lambda2 = 0.2)
wls_l1_nuc.cov <- function(y, x, lambda1, lambda2, weights = NULL, thresh = 1e-5,
                       maxit = 100, mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                       trace.it = F, offset = F, max.rank = 20) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  q <- ncol(x)
  y0 <- y
  omega <- !is.na(y)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) alpha0 <- rep(0, q)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  theta <- theta0
  objective <- NULL
  y <- y0
  y[!omega] <- 0
  error <- 1
  iter <- 0
  if(is.null(weights)) weights <- matrix(rep(1, n*p), nrow = n)
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    y.tmp <- mu + alpha.mat + theta
    y <- y0
    y[!omega] <- y.tmp[!omega]
    mu.tmp <- mu
    alpha.tmp <- alpha
    theta.tmp <- theta
    if(offset) mu <-  sum(weights * (y - alpha.mat - theta.tmp), na.rm = T) / sum(weights) else mu <- 0
    alpha <- glmnet(x, c(y - mu - theta.tmp), family = "gaussian", lambda = seq(1e3*lambda2, lambda2, length.out = 100),
                    intercept = FALSE, thresh = thresh, weights = c(weights))$beta[, 100]
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    svd_theta <- wlra(y - mu - alpha.mat, w = weights, lambda = lambda1, x0 = NULL,
                              thresh = thresh, rank.max = max.rank)
    u <- svd_theta$u
    d <- svd_theta$d
    v <- svd_theta$v
    if(is.null(dim(u))){
      theta <- d * u%*%t(v)
    } else {
      theta <- u%*%diag(d)%*%t(v)
    }
    objective <- c(objective, sum((1 / 2) * weights * (y0 - mu - alpha.mat - theta)^2, na.rm = T) + lambda1 * sum(d) +
                     lambda2 * sum(abs(alpha)))
    if(iter == 1) {
      error <- 1
    } else{
      error <- abs(objective[iter]-objective[iter - 1])/abs(objective[iter])
    }
    if(trace.it && (iter%%10 == 0)){
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  y <- mu + alpha.mat + theta
  y.imputed <- y0
  y.imputed[is.na(y0)] <- y[is.na(y0)]
  return(list(y = y0, y.imputed = y.imputed, param = mu + alpha.mat + theta, mu = mu,
              alpha = alpha, theta = theta, objective = objective, iter = iter))

}

#' rwls_l1_nuc.cov
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
#' \item{y}{the original data matrix}
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' @export
#' @import FactoMineR
#'
#' @examples
#' n = 6; p = 2
#' y1 <- matrix(rnorm(mean = 0, n * p), nrow = n)
#' y2 <- matrix(rnorm(mean = 0, n * p), nrow = n)
#' y3 <- matrix(rnorm(mean = 2, n * p), nrow = n)
#' y <- cbind(matrix(rnorm(mean = c(y1), n * p), nrow = n),
#'            matrix(rbinom(n * p, prob = c(exp(y2)/(1+exp(y2))), size = 1), nrow = n),
#'            matrix(rpois(n * p, lambda = c(exp(y3))), nrow = n))
#' var.type <- c(rep("gaussian", p), rep("binary", p), rep("poisson", p))
#' idx_NA <- sample(1:(3 * n * p), size = round(0.1 * 3 * n * p))
#' y[idx_NA] <- NA
#' x <- matrix(rnorm(6*6*2), nrow = 6*6)
#' nl <- rep(1, 6)
#' res <- rwls_l1_nuc.cov(y, x, var.type, 0.1, 0.2, nl)
rwls_l1_nuc.cov <- function(y, x, var.type, lambda1, lambda2, nlevel = NULL, maxit = 100,
                            mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                            thresh = 1e-5, trace.it = F, offset = F, scale = F, max.rank = 20,
                            vt2 = NULL, wmax = NULL){
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  q <- ncol(x)
  x <- as.matrix(x)
  x <- matrix(as.numeric(x), nrow = nrow(x))
  omega <- !is.na(y)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) alpha0 <- rep(0, q)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
  theta <- theta0
  param <- mu + alpha.mat + theta
  objective <- NULL
  y0 <- y
  error <- 1
  iter <- 0
  if(scale == T){
    moy <- colMeans(y, na.rm = T)
    moy[vt2 == "binary"] <- log(moy[vt2 == "binary"]/(1-moy[vt2 == "binary"]))
    moy[vt2 == "poisson"] <- log(moy[vt2 == "poisson"])
    moy[vt2 == "categorical"] <- log(moy[vt2 == "categorical"]/(1-moy[vt2 == "categorical"]))
    scaling <- rep(0, p)
    p1 <- sum(var.type == "gaussian")
    p2 <- sum(var.type == "binary")
    p3 <- sum(var.type == "poisson")
    p4 <- sum(var.type == "categorical")
    if(p1>0){
      scaling[vt2 == "gaussian"] <- 0.5*colSums((y[, vt2 == "gaussian"] - t(matrix(rep(moy[vt2 == "gaussian"],
                                                                                       n), nrow = p1)))^2,na.rm = T)/(n-1)
    }
    if(p2>0){
      scaling[vt2 == "binary"] <- colSums(- (y[, vt2 == "binary"] * t(matrix(rep(moy[vt2 == "binary"],
                                                                                 n), nrow = p2))) +
                                            log(1 + exp(t(matrix(rep(moy[vt2 == "binary"],
                                                                     n), nrow = p2)))),na.rm = T)/(n-1)

    }
    if(p3>0){
      scaling[vt2 == "poisson"] <- colSums(- (y[, vt2 == "poisson"] * t(matrix(rep(moy[vt2 == "poisson"],
                                                                                   n), nrow = p3))) +
                                             exp(t(matrix(rep(moy[vt2 == "poisson"],n), nrow = p3))) + log_factorial(y[, vt2 == "poisson"]),na.rm = T)/(n-1)

    }
    if(p4>0){
      scaling[vt2 == "categorical"] <- colSums(- (y[, vt2 == "categorical"] * t(matrix(rep(moy[vt2 == "categorical"],
                                                                                           n), nrow = sum(vt2 == "categorical")))) +
                                                 log(1 + exp(t(matrix(rep(moy[vt2 == "categorical"],
                                                                          n), nrow = sum(vt2 == "categorical"))))),na.rm = T)/(n-1)

    }
  } else scaling <- rep(1, ncol(y))
  sc <- matrix(rep(scaling, n), nrow = n, byrow = T)
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    mu.tmp <- mu
    alpha.tmp <- alpha
    theta.tmp <- theta
    yv <- quad_approx(y0, param, var.type, nlevel, wmax)
    ytilde <- yv$ytilde
    vtilde2 <- yv$vtilde2
    lambda1w <- lambda1 / (max(vtilde2))
    lambda2w <- lambda2 / (max(vtilde2))
    vtilde2 <- vtilde2 / max(vtilde2)
    if(scale==T){
      vtilde2 <- sweep(vtilde2, 2, scaling, "/")
    }
    res_approx <- wls_l1_nuc.cov(ytilde, x, lambda1w, lambda2w, weights = vtilde2,
                                 thresh = thresh, mu0 = mu.tmp, alpha0 = alpha.tmp,
                                 theta0 = theta.tmp, trace.it = F, maxit = maxit,
                                 offset = offset)
    mu <- res_approx$mu
    alpha <- res_approx$alpha
    theta <- res_approx$theta
    alpha.mat <- matrix(matrix(as.numeric(x), nrow=n*p)%*%alpha, nrow=n)
    param <- mu + alpha.mat + theta
    res <- bls.cov(y0, x, mu, alpha, theta, mu.tmp, alpha.tmp, theta.tmp,
                   b=0.5, lambda1, lambda2, var.type, thresh, sc=sc)
    mu <- res$mu
    alpha <- res$alpha
    theta <- res$theta
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    param <- mu + alpha.mat + theta
    gaus <- (1/2) * sum(sc[, var.type == "gaussian"]*(y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum(sc[, var.type == "poisson"]*(- (y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                                               exp(param[, var.type == "poisson"])), na.rm = T)
    binom <- sum(sc[, var.type == "binary"]*(- (y0[, var.type == "binary"] * param[, var.type == "binary"]) +
                                               log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if(sum(var.type=="categorical")>0){
      for(j in 1:sum(var.type=="categorical")){
        tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
        truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum(sc[, vt2 == "categorical"]*(- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                  log(truc))), na.rm = T)
    } else cat <- 0

    d <- svd(theta)$d
    if(length(lambda2) == 1){
      objective <- c(objective, min(.Machine$double.xmax,
                                    (pois + gaus + binom + cat + lambda1 * sum(d) + lambda2 * sum(abs(alpha)))))
    } else {
      objective <- c(objective, min(.Machine$double.xmax,
                                    (pois + gaus + binom + cat + lambda1 * sum(d) + sum(lambda2 * t(abs(alpha))))))
    }

    if(iter == 1) {
      error <- 1
    } else{
      if(objective[iter] == .Machine$double.xmax) {
        error <- 1
      } else error <- abs(objective[iter]-objective[iter - 1])/abs(objective[iter])

    }
    if(all(var.type == "gaussian")) error <- 0
    if(trace.it){
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  if(error < thresh) cvg = T else cvg = F
  y <- mu + alpha.mat + theta
  y[, var.type == "poisson"] <- exp(y[, var.type == "poisson"])
  y[, var.type == "binary"] <- exp(y[, var.type == "binary"])/(1+exp(y[, var.type == "binary"]))
  y.imputed <- y0
  y.imputed[is.na(y0)] <- y[is.na(y0)]
  return(list(y = y0, y.imputed = y.imputed, param = mu + alpha.mat + theta, mu = mu,
              alpha = alpha, theta = theta, objective = objective))
}


