#' rwls_l1_nuc.lr
#'
#' @param y nxp observation matrix
#' @param var.type vector of size p indicating types of columns in y (gaussian, binary, poisson)
#' @param lambda1 positive number, value of the nuclear norm regularization parameter
#' @param nlevel vector of integers indicating the number of levels of each factor in y
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param upper real number, upper bound on entries of theta
#' @param lower real number, lower bound on entries of theta
#' @param theta0 matrix of size (nb of ind.)x(number of variables), initial value of the individual effect, default 0
#' @param trace.it boolean, if TRUE information about convergence will be displayes, default FALSE
#' @param offset boolean, if TRUE offset is computed, otherwise set to 0, default FALSE
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param max.rank integer, maximum rank of interaction matrix theta
#' @param vt2 vector indicating types of the columns of the extended data frame y (with dummies for every category)

#' @return A list vith the following elements
#' \item{y}{the original data matrix}
#' \item{param}{a (nb individuals) x (nb variables) parameter matrix}
#' \item{theta}{a (nb individuals) x (nb variables) parameter matrix}
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
#' nl <- rep(1, 6)
#' res <- rwls_l1_nuc.lr(y, var.type, 0.1, nl)
rwls_l1_nuc.lr <- function(y, var.type, lambda1, nlevel = NULL, maxit = 100, upper = 12,
                            lower = -12, theta0 = NULL, thresh = 1e-5, trace.it = F,
                            offset = F, scale = F, max.rank = 20, vt2 = NULL, wmax = NULL){
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  omega <- !is.na(y)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  theta <- theta0
  param <- theta
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
    theta.tmp <- theta
    param.tmp <- param
    yv <- quad_approx(y0, param, var.type, nlevel, wmax)
    ytilde <- yv$ytilde
    vtilde2 <- yv$vtilde2
    lambda1w <- lambda1 / (max(vtilde2))
    vtilde2 <- vtilde2 / max(vtilde2)
    if(scale==T){
      vtilde2 <- sweep(vtilde2, 2, scaling, "/")
    }
    svd_theta <- wlra(x = ytilde, w = vtilde2, lambda = lambda1w, x0 = NULL, thresh = thresh,
                      rank.max = max.rank)
    u <- svd_theta$u
    d <- svd_theta$d
    v <- svd_theta$v
    if(is.null(dim(u))){
      theta <- d * u%*%t(v)
    } else {
      theta <- u%*%diag(d)%*%t(v)
    }
    theta[theta > upper] <- upper
    theta[theta < lower] <- lower
    param <- theta
    res <- bls.lr(y0, theta, theta.tmp, b = 0.5, lambda1, var.type, thresh)
    theta <- res$theta
    param <- theta
    gaus <- (1 / 2) * sum(sc[, var.type == "gaussian"]*(y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
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
    objective <- c(objective, min(.Machine$double.xmax, (pois + gaus + binom + cat + lambda1 * sum(d))))

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
  y <- theta
  y[, var.type == "poisson"] <- exp(y[, var.type == "poisson"])
  y[, var.type == "binary"] <- exp(y[, var.type == "binary"])/(1+exp(y[, var.type == "binary"]))
  y.imputed <- y0
  y.imputed[is.na(y.imputed)] <- y[is.na(y.imputed)]
  return(list(y = y0, theta = theta, objective = objective, y.imputed = y.imputed))
}


