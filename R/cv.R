#' cv.mimi
#' Perform cross-validation for mimi model with covariates/groups along regularization path
#' @param y nxp matrix of observations
#' @param model either one of "groups" or "covariates", indicating which model should be fitted
#' @param x (np)xN matrix of covariates (optional)
#' @param groups factor of length n indicating groups
#' @param nb.boot integer indicating the number of bootstrap samples for cross-validation
#' @param prob number in (0,1) proportion of missing entries added in the bootstrap samples
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param var.type vector of length p indicating the types of the columns of y (gaussian, binomial, poisson)
#' @param lambda1.max positive number, max value of nuclear norm regularization parameter
#' @param lambda2.max positive number, max value of l1 norm regularization parameter
#' @param lambda1.min positive number, min value of nuclear norm regularization parameter
#' @param lambda2.min positive number, min value of l1 norm regularization parameter
#' @param size.grid integer, length of the grid of lambda1 and lambda2 for cross-validation
#' @param length integer, length of the grid of lambda1 and lambda2 for warm start
#' @param mu0 real number, optional starting point for offset
#' @param alpha0 Nxp matrix, optional starting point for group effects
#' @param theta0 nxp matrix, optional starting point for interaction matrix
#' @param trace.it boolean, whether convergence information should be printed
#' @param offset boolean, whether an offset should be fitted
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 used for warm start (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 used for warm start (log scale)}
#' \item{list.res}{list of results of rwls_l1_nuc (same length as grid of lambda)}
#' \item{list.mu}{list of offsets (same length as grid of lambda)}
#' \item{list.alpha}{list of regression parameters (same length as grid of lambda)}
#' \item{list.theta}{list of interaction matrices (same length as grid of lambda)}
#' @export
#' @import softImpute
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
#' res <- cv.mimi(y, model = "cov", x = x, var.type = var.type, maxit = 2, nb.boot = 2, size.grid = 2)
cv.mimi <- function(y, model = c("groups", "covariates"), x = NULL, groups = NULL,
                    var.type, lambda1.min = NULL, lambda2.min = NULL,
                    lambda1.max = NULL, lambda2.max = NULL, maxit = 100,
                    mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                    thresh = 1e-4, trace.it = F, length = 20, offset = F, nb.boot = 10, prob = 0.1,
                    size.grid = 20)
{
  if(model == "groups"){
    return (cv.mimi.multi(y, groups, nb.boot, prob , thresh,
                          maxit, var.type, lambda1.max, lambda2.max,
                          lambda1.min, lambda2.min, size.grid, length,
                          mu0, alpha0, theta0, trace.it, offset))
  } else{
    return(cv.mimi.cov(y, x, nb.boot, prob, thresh,
                       maxit, var.type, lambda1.max,
                       lambda2.max, lambda1.min,
                       lambda2.min, size.grid, length, mu0,
                       alpha0, theta0, trace.it, offset))
  }
}




#' cv.mimi.multi
#'
#' @param y nxp observations matrix
#' @param groups factor of length n indicating groups
#' @param nb.boot integer indicating the number of bootstrap samples for cross-validation
#' @param prob number in (0,1) proportion of missing entries added in the bootstrap samples
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param var.type vector of length p indicating the types of the columns of y (gaussian, binomial, poisson)
#' @param lambda1.max positive number, max value of nuclear norm regularization parameter
#' @param lambda2.max positive number, max value of l1 norm regularization parameter
#' @param lambda1.min positive number, min value of nuclear norm regularization parameter
#' @param lambda2.min positive number, min value of l1 norm regularization parameter
#' @param size.grid integer, length of the grid of lambda1 and lambda2 for cross-validation
#' @param length integer, length of the grid of lambda1 and lambda2 for warm start
#' @param mu0 real number, optional starting point for offset
#' @param alpha0 Nxp matrix, optional starting point for group effects
#' @param theta0 nxp matrix, optional starting point for interaction matrix
#' @param trace.it boolean, whether convergence information should be printed
#' @param offset boolean, whether an offset should be fitted

#'
#' @return A list containing the following elements
#' \item{lambda1}{The selected value of lambda1 (regularization parameter for nuclear norm)}
#' \item{lambda2}{The selected value of lambda2 (regularization parameter for l1 norm)}
#' \item{estim.cv}{An object resulting from the mimi method, estimated with the selected values of lambda1 and lambda2}
#' \item{lambda1.grid}{A vector of size length containing the grid of lambda1 used for cv}
#' \item{lambda2.grid}{A vector of size length containing the grid of lambda2 used for cv}
#' @export
#' @import softImpute
#'
#' @examples
#' y <- matrix(rnorm(4 * 3), nrow = 4)
#' y[sample(1:12, size = 3)] <- NA
#' groups <- c(1,1,2,2)
#' groups <- as.factor(groups)
#' vt <- rep("gaussian", 3)
#' res <- cv.mimi.multi(y, groups, 2, var.type = vt, maxit = 1, size.grid = 2)
cv.mimi.multi <- function(y, groups, nb.boot = 10, prob = 0.1, thresh = 1e-4,
                          maxit = 100, var.type, lambda1.max = NULL, lambda2.max = NULL,
                          lambda1.min = NULL, lambda2.min = NULL, size.grid = 20, length = 20,
                          mu0 = NULL, alpha0 = NULL, theta0 = NULL, trace.it = F,
                          offset = F){
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y)
  w.groups <- aggregate(matrix(omega, nrow = n), list(groups), sum)[, 2:(p + 1)]
  w.groups[w.groups == 0] <- 1
  y.mean <-  sum(y, na.rm = T) / sum(omega)
  y.center <- y.center <- aggregate(y - y.mean, list(groups), sum, na.rm = T)[, 2:(p+1)] / w.groups
  if(is.null(lambda2.max)) lambda2.max <- 2 * max(w.groups * y.center)
  ncenters <- aggregate(rep(1, n), list(groups), sum)[,2]
  N <- nlevels(groups)
  warm.start <- softImpute(y - y.mean - matrix(rep(as.matrix(y.center), rep(ncenters, p)), nrow = n))
  if(is.null(lambda1.max)) lambda1.max <- 10*max(warm.start$d)
  if(is.null(lambda1.min)) lambda1.min <- lambda1.max / 1e3
  if(is.null(lambda2.min)) lambda2.min <- lambda2.max / 1e3
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1.min), length.out = size.grid)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2.min), length.out = size.grid)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) {
    alpha0 <- rep(0, N * p)
    alpha.rep <- matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
    alpha0 <- matrix(as.matrix(alpha0), nrow = N)
  } else{
    alpha.rep <- matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
    alpha <- matrix(as.matrix(alpha0), nrow = N)
  }
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  alpha.rep <- matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
  theta <- theta0
  iter <- 1
  obs.idx <- which(omega)
  na.func <- function(x, prob){
    yp <- y
    yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    yp
  }
  y.list <- lapply(1:nb.boot, function(i) na.func(y, prob))
  pred.list <- lapply(y.list, function(yy) (is.na(yy) * (!is.na(y))))
  res.cv <- list()
  for(k in 1:nb.boot){
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) mimi.multi(y.list[[k]], groups = groups, var.type = var.type,
                                                                                                lambda1 = exp(lambda1.grid.log[i]),
                                                                                                lambda2 = exp(lambda2.grid.log[j]),length = length,
                                                                                                thresh = thresh, trace.it = trace.it, maxit = maxit)))
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) res.cv[[k]][[i]][[j]]$list.res[[length]]))

    res.cv[[k]] <- sapply(1:length(lambda1.grid.log), function(i) sapply(1:length(lambda2.grid.log),
                                                                         function(j) sum((res.cv[[k]][[i]][[j]]$y.imputed * pred.list[[k]] - y * pred.list[[k]])^2, na.rm = T)))


  }

  error <- Reduce('+', res.cv)/nb.boot
  idx <- which(error == min(error, na.rm =T), arr.ind = TRUE)
  lambda1.cv <- exp(lambda1.grid.log)[idx[1]]
  lambda2.cv <- exp(lambda2.grid.log)[idx[2]]
  estim.cv <- mimi.multi(y, groups = groups, lambda1 = lambda1.cv, lambda2 = lambda2.cv,
                         var.type = var.type, trace.it = trace.it, length = length,
                         thresh = thresh, maxit = maxit)
  return(list(lambda1 = exp(lambda1.cv), lambda2 = exp(lambda2.cv), estim.cv = estim.cv,
              lambda1.grid = exp(lambda1.grid.log), lambda2.grid = exp(lambda2.grid.log)))

}


#' cv.cov
#'
#' @param y nxp observation matrix
#' @param x (np)xq matrix of covariates
#' @param nb.boot integer indicating the number of bootstrap samples for cross-validation
#' @param prob number in (0,1) proportion of missing entries added in the bootstrap samples
#' @param thresh positive number, convergence criterion
#' @param maxit integer, maximum number of iterations
#' @param var.type vector of length p indicating the types of the columns of y (gaussian, binomial, poisson)
#' @param lambda1.max positive number, max value of nuclear norm regularization parameter
#' @param lambda2.max positive number, max value of l1 norm regularization parameter
#' @param lambda1.min positive number, min value of nuclear norm regularization parameter
#' @param lambda2.min positive number, min value of l1 norm regularization parameter
#' @param size.grid integer, length of the grid of lambda1 and lambda2 for cross-validation
#' @param length integer, length of the grid of lambda1 and lambda2
#' @param mu0 real number, optional starting point for offset
#' @param alpha0 Nxp matrix, optional starting point for group effects
#' @param theta0 nxp matrix, optional starting point for interaction matrix
#' @param trace.it boolean, whether convergence information should be printed
#' @param offset boolean, whether an offset should be fitted
#'
#' @return A list containing the following elements
#' \item{lambda1}{The selected value of lambda1 (regularization parameter for nuclear norm)}
#' \item{lambda2}{The selected value of lambda2 (regularization parameter for l1 norm)}
#' \item{estim.cv}{An object resulting from the mimi method, estimated with the selected values of lambda1 and lambda2}
#' \item{lambda1.grid}{A vector of size length containing the grid of lambda1 used for cv}
#' \item{lambda2.grid}{A vector of size length containing the grid of lambda2 used for cv}
#' @export
#' @import softImpute
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
#' res <- cv.mimi.cov(y, x, var.type = var.type, maxit = 2, nb.boot = 2, thresh =1)
cv.mimi.cov <- function(y, x, nb.boot = 5, prob = 0.1, thresh = 1e-4,
                        maxit = 100, var.type, lambda1.max = NULL,
                        lambda2.max = NULL, lambda1.min = NULL,
                        lambda2.min = NULL, size.grid = 10, length = 20, mu0 = NULL,
                        alpha0 = NULL, theta0 = NULL, trace.it = F, offset = T){
  x <- as.matrix(x)
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  q <- ncol(x)
  omega <- !is.na(y)
  if(is.null(lambda2.max)) lambda2.max <- 2 * max(y, na.rm = T)
  warm.start <- softImpute(as.matrix(y))
  if(is.null(lambda1.max)) lambda1.max <- max(warm.start$d)
  if(is.null(lambda1.min)) lambda1.min <- lambda1.max / (3 * 1e3)
  if(is.null(lambda2.min)) lambda2.min <- lambda2.max / (3 * 1e3)
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1.min), length.out = size.grid)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2.min), length.out = size.grid)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) alpha0 <- rep(0, q)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
  theta <- theta0
  iter <- 1
  obs.idx <- which(omega)
  if(is.null(nb.boot)) nb.boot <- 10
  na.func <- function(y, prob){
    yp <- y
    yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    yp
  }
  y.list <- lapply(1:nb.boot, function(i) na.func(y, prob))
  pred.list <- lapply(y.list, function(yy) (is.na(yy) * (!is.na(y))))
  res.cv <- list()
  for(k in 1:nb.boot){
    res.cv[[k]] <- lapply(1:size.grid,
                          function(i) lapply(1:size.grid,
                                             function(j) mimi.cov(y.list[[k]], x,
                                                                   var.type = var.type,
                                                                   lambda1 = exp(lambda1.grid.log[i]),
                                                                   lambda2 = exp(lambda2.grid.log[j]),
                                                                   thresh = thresh, trace.it = trace.it,
                                                                   offset = offset, lambda1.max = lambda1.max,
                                                                   lambda2.max = lambda2.max)))

    res.cv[[k]] <- lapply(1:size.grid,
                          function(i) lapply(1:size.grid,
                                             function(j) res.cv[[k]][[i]][[j]]$list.res[[length]]))


    res.cv[[k]] <- sapply(1:size.grid,
                          function(i) sapply(1:size.grid,
                                             function(j) sum((res.cv[[k]][[i]][[j]]$y.imputed * pred.list[[k]] - y * pred.list[[k]])^2, na.rm = T)))


  }

  error <- Reduce('+', res.cv)/nb.boot
  idx <- which(error == min(error, na.rm =T), arr.ind = TRUE)
  lambda1.cv <- exp(lambda1.grid.log)[idx[1]]
  lambda2.cv <- exp(lambda2.grid.log)[idx[2]]
  estim.cv <- mimi.cov(y, x, lambda1 = lambda1.cv, lambda2 = lambda2.cv,
                        var.type = var.type, trace.it = T, length = length, thresh = thresh)
  return(list(lambda1 = lambda1.cv, lambda2 = lambda2.cv, estim.cv = estim.cv, error = error,
              lambda1.grid = exp(lambda1.grid.log), lambda2.grid = exp(lambda2.grid.log)))

}
