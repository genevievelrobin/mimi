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
                       trace.it = F, offset = F) {
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
    alpha <- glmnet(x, c(y), family = "gaussian", lambda = seq(1e3*lambda2, lambda2, length.out = 100),
                    intercept = FALSE, thresh = thresh, weights = c(weights))$beta[, 100]
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    svd_theta <- wlra(y - mu - alpha.mat, w = weights, lambda = lambda1, x0 = NULL,
                              thresh = thresh)
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
#' @param var.type vector of size p indicating types of columns in y (gaussian, binomial, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param lambda2 positive number, value of the l1 norm regularization parameter
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param upper real number, upper bound on entries of mu, alpha and theta
#' @param lower real number, lower bound on entries of mu, alpha and theta
#' @param mu0 real number, initial value of the offset, default 0
#' @param alpha0 matrix of size (nb of groups)x(number of variables), initial value of the group effect, default 0
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param scale boolean, indicates whether or not cost functions should be scaled by column
#' @param offset boolean, if TRUE offset is computed, otherwise set to 0, default FALSE
#'
#' @return A list vith the following elements
#' \item{y}{the original data matrix}
#' \item{y.imputed}{the original data matrix where missing entries are imputed by their estimated means}
#' \item{param}{the estimated parameter matrix}
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta}{a (nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' @export
#'
#' @examples
#' n = 6; p = 2
#' y1 <- matrix(rnorm(mean = 0, n * p), nrow = n)
#' y2 <- matrix(rnorm(mean = 0, n * p), nrow = n)
#' y3 <- matrix(rnorm(mean = 2, n * p), nrow = n)
#' y <- cbind(matrix(rnorm(mean = c(y1), n * p), nrow = n),
#'            matrix(rbinom(n * p, prob = c(exp(y2)/(1+exp(y2))), size = 1), nrow = n),
#'            matrix(rpois(n * p, lambda = c(exp(y3))), nrow = n))
#' var.type <- c(rep("gaussian", p), rep("binomial", p), rep("poisson", p))
#' idx_NA <- sample(1:(3 * n * p), size = round(0.1 * 3 * n * p))
#' y[idx_NA] <- NA
#' x <- matrix(rnorm(6*6*2), nrow = 6*6)
#' res <- rwls_l1_nuc.cov(y, x, var.type = var.type, lambda1 = 0.1, lambda2 = 0.2)
rwls_l1_nuc.cov <- function(y, x, var.type, lambda1, lambda2, maxit = 100, upper = 12,
                        lower = -12, mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                        thresh = 1e-5, trace.it = F, scale = F, offset = F){
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  q <- ncol(x)
  if(scale == T){
    moy <- colMeans(y, na.rm = T)
    moy[var.type == "binomial"] <- log(moy[var.type == "binomial"]/(1-moy[var.type == "binomial"]))
    moy[var.type == "poisson"] <- log(moy[var.type == "poisson"])
    scaling <- rep(0, p)
    p1 <- sum(var.type == "gaussian")
    p2 <- sum(var.type == "binomial")
    p3 <- sum(var.type == "poisson")
    if(p1>0){
      scaling[var.type == "gaussian"] <- 0.5*colSums((y[, var.type == "gaussian"] - t(matrix(rep(moy[var.type == "gaussian"],
                                                                                                 n), nrow = p1)))^2,na.rm = T)/(n-1)
    }
    if(p2>0){
      scaling[var.type == "binomial"] <- colSums(- (y[, var.type == "binomial"] * t(matrix(rep(moy[var.type == "binomial"],
                                                                                               n), nrow = p2))) +
                                                   log(1 + exp(t(matrix(rep(moy[var.type == "binomial"],
                                                                            n), nrow = p2)))),na.rm = T)/(n-1)

    }
    if(p3>0){
      scaling[var.type == "poisson"] <- colSums(- (y[, var.type == "poisson"] * t(matrix(rep(moy[var.type == "poisson"],
                                                                                             n), nrow = p3))) +
                                                  exp(t(matrix(rep(moy[var.type == "poisson"],n), nrow = p3))) + log_factorial(y[, var.type == "poisson"]),na.rm = T)/(n-1)

    }
  }
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
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    mu.tmp <- mu
    alpha.tmp <- alpha
    theta.tmp <- theta
    yv <- lapply(1:p, function(j) quad_approx(y0[, j], param[, j], var.type[j]))
    ytilde <- do.call(cbind, lapply(yv, function(t) t$ytilde))
    vtilde2 <- do.call(cbind, lapply(yv, function(t) t$vtilde2))
    if(scale==T){
      vtilde2 <- sweep(vtilde2, 2, scaling, "/")
    }
    lambda1w <- lambda1 / (max(vtilde2))
    lambda2w <- lambda2 / (max(vtilde2))
    vtilde2 <- vtilde2 / max(vtilde2)
    res_approx <- wls_l1_nuc.cov(ytilde, x, lambda1w, lambda2w, weights = vtilde2,
                                 thresh = thresh, mu0 = mu.tmp, alpha0 = alpha.tmp,
                                 theta0 = theta.tmp, trace.it = F, maxit = maxit,
                                 offset = offset)
    mu <- res_approx$mu
    alpha <- res_approx$alpha
    theta <- res_approx$theta
    if(mu > upper) mu <- upper
    if(mu < lower) mu <- lower
    alpha[alpha > upper] <- upper
    theta[theta > upper] <- upper
    alpha[alpha < lower] <- lower
    theta[theta < lower] <- lower
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    param <- mu + alpha.mat + theta
    res <- bls.cov(y0, x, mu, alpha, theta, mu.tmp, alpha.tmp, theta.tmp,
               b = 0.5, lambda1, lambda2, var.type, thresh)
    mu <- res$mu
    alpha <- res$alpha
    theta <- res$theta
    alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
    param <- mu + alpha.mat + theta
    gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum(- (y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                  exp(param[, var.type == "poisson"]), na.rm = T)
    binom <- sum(- (y0[, var.type == "binomial"] * param[, var.type == "binomial"]) +
                   log(1 + exp(param[, var.type == "binomial"])), na.rm = T)
    d <- svd(theta)$d
    if(length(lambda2) == 1){
      objective <- c(objective, min(.Machine$double.xmax,
                                    (pois + gaus + binom + lambda1 * sum(d) + lambda2 * sum(abs(alpha)))))
    } else {
      objective <- c(objective, min(.Machine$double.xmax,
                                    (pois + gaus + binom + lambda1 * sum(d) + sum(lambda2 * t(abs(alpha))))))
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
  y[, var.type == "binomial"] <- exp(y[, var.type == "binomial"])/(1+exp(y[, var.type == "binomial"]))
  y.imputed <- y0
  y.imputed[is.na(y0)] <- y[is.na(y0)]
  return(list(y = y0, y.imputed = y.imputed, param = mu + alpha.mat + theta, mu = mu,
              alpha = alpha, theta = theta, objective = objective))
}


