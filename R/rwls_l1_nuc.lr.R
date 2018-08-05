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
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    theta.tmp <- theta
    param.tmp <- param
    yv <- quad_approx(y0, param, var.type, nlevel, wmax)
    ytilde <- yv$ytilde
    vtilde2 <- yv$vtilde2
    if(scale==T){
      scaling <- apply(ytilde, 2, sd, na.rm = T)
      ytilde <- sweep(ytilde, 2, scaling, "/")
    }
    lambda1w <- lambda1 / max(vtilde2)
    vtilde2 <- vtilde2 / max(vtilde2)
    svd_theta <- wlra(x = ytilde, w = vtilde2, lambda = lambda1w, x0 = NULL, thresh = 0.1*thresh,
                      rank.max = max.rank, maxit = 100)
    u <- svd_theta$u
    d <- svd_theta$d
    v <- svd_theta$v
    if(is.null(dim(u))){
      theta <- d * u%*%t(v)
    } else {
      theta <- u%*%diag(d)%*%t(v)
    }
    if(sum(var.type == "categorical")) theta[, vt2 == "categorical"] <- sweep(theta[, vt2 == "categorical"], 1, rowMeans(theta[, vt2 == "categorical"]))
    theta[theta > upper] <- upper
    theta[theta < lower] <- lower
    if(scale) theta <- sweep(theta, 2, scaling, "*")
    res <- bls.lr(y0, theta, theta.tmp, b = 0.5, lambda1, vt2, thresh)
    theta <- res$theta
    param <- theta
    #obj <- sum(vtilde2*(ytilde-sweep(theta, 2, scaling, "/"))^2, na.rm =T) + lambda1w*sum(svd(sweep(theta, 2, scaling, "/"))$d)
    gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - theta[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum((- (y0[, var.type == "poisson"] * theta[, var.type == "poisson"]) +
                   exp(theta[, var.type == "poisson"])), na.rm = T)
    binom <- sum((- (y0[, var.type == "binary"] * theta[, var.type == "binary"]) +
                    log(1 + exp(theta[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if(sum(var.type=="categorical")>0){
      for(j in 1:sum(var.type=="categorical")){
        tt <- rowSums(exp(theta[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
        truc <- cbind(truc, matrix(rep(theta, nlevel[which(var.type=="categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum((- (y0[, vt2 == "categorical"] * (theta[, vt2 == "categorical"]) -
                       log(truc))), na.rm = T)
    } else cat <- 0
    d <- svd(theta)$d
    objective <- c(objective, min(.Machine$double.xmax, (pois + gaus + binom + cat + lambda1 * sum(d))))
    #objective <- c(objective, min(.Machine$double.xmax, obj))
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
  y[, vt2 == "poisson"] <- exp(y[, vt2 == "poisson"])
  y[, vt2 == "binary"] <- exp(y[, vt2 == "binary"])/(1+exp(y[, vt2 == "binary"]))
  if(sum(var.type == "categorical")>0){
    for(j in 1:sum(var.type == "categorical")){
      count <- sum(nlevel[1:(which(var.type == "categorical")[j]-1)])
      y[, (count+1):(count+nlevel[which(vt2 == "categorical")[j]])] <- t(sapply(1:n, function(i) exp(y[i, (count+1):(count+nlevel[which(vt2 == "categorical")[j]])])/sum(exp(y[i, (count+1):(count+nlevel[which(vt2 == "categorical")[j]])]))))
    }
  }
  y.imputed <- y0
  y.imputed[is.na(y.imputed)] <- y[is.na(y.imputed)]
  return(list(y = y0, theta = theta, objective = objective, y.imputed = y.imputed))
}


