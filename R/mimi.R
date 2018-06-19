#' mimi
#' Compute solution of mimi with covariates/groups along regularization path
#' @param y nxp matrix of observations
#' @param model either one of "groups" or "covariates", indicating which model should be fitted
#' @param x (np)xN matrix of covariates (optional)
#' @param groups factor of length n indicating groups (optional)
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binomial or poisson)
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param maxit integer maximum number of iterations
#' @param mu0 real number initial value of offset (optional)
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param lambda1.max positive number starting value of regularization parameter
#' @param lambda2.max positive number starting value of regularization parameter
#' @param length length of the grid of regularization parameters
#' @param upper real number upper bound on entries of alpha and theta
#' @param lower real number lower bound on entries of alpha and theta
#' @param scale boolean indicating whether or not the cost functions should be scaled by column
#' @param offset boolean indicating whether or not an offset should be fitted
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 used for warm start (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 used for warm start (log scale)}
#' \item{list.res}{list of results of rwls_l1_nuc (same length as grid of lambda)}
#' \item{list.mu}{list of offsets (same length as grid of lambda)}
#' \item{list.alpha}{list of regression parameters (same length as grid of lambda)}
#' \item{list.theta}{list of interaction matrices (same length as grid of lambda)}
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
#' res <- mimi(y, model = "covariates", x = x, var.type = var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi <- function(y, model = c("groups", "covariates"), x = NULL, groups = NULL,
                 var.type, lambda1, lambda2, maxit = 100,
                 mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                 thresh = 1e-4, trace.it = F, lambda1.max = NULL,
                 lambda2.max = NULL, length = 20, upper = 12,
                 lower = -12, scale = F, offset = F)
{
  if(model == "groups"){
    return (mimi.multi(y, groups, var.type, lambda1, lambda2, maxit,
                       mu0, alpha0, theta0, thresh, trace.it,
                       lambda1.max, lambda2.max, length, upper,
                       lower, scale, offset))
  } else{
    return(mimi.cov(y, x, var.type, lambda1, lambda2, maxit,
                      mu0, alpha0, theta0, thresh, trace.it,
                      lambda1.max, lambda2.max, length, upper,
                      lower, scale, offset))
  }
}



#' mimi.cov
#' Compute solution of mimi with covariates effects along regularization path
#' @param y nxp matrix of observations
#' @param x (np)xN matrix of covariates
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binomial or poisson)
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param maxit integer maximum number of iterations
#' @param mu0 real number initial value of offset (optional)
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param lambda1.max positive number starting value of regularization parameter
#' @param lambda2.max positive number starting value of regularization parameter
#' @param length length of the grid of regularization parameters
#' @param upper real number upper bound on entries of alpha and theta
#' @param lower real number lower bound on entries of alpha and theta
#' @param scale boolean indicating whether or not the cost functions should be scaled by column
#' @param offset boolean indicating whether or not an offset should be fitted
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 used for warm start (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 used for warm start (log scale)}
#' \item{list.res}{list of results of rwls_l1_nuc (same length as grid of lambda)}
#' \item{list.mu}{list of offsets (same length as grid of lambda)}
#' \item{list.alpha}{list of regression parameters (same length as grid of lambda)}
#' \item{list.theta}{list of interaction matrices (same length as grid of lambda)}
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
#' res <- mimi.cov(y, x, var.type = var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi.cov <- function(y, x, var.type, lambda1, lambda2, maxit = 1e3,
                     mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                     thresh = 1e-4, trace.it = F, lambda1.max = NULL,
                     lambda2.max = NULL, length = 20, upper = 12,
                     lower = -12, scale = F, offset = F)
{
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y)
  q <- ncol(x)
  if(is.null(lambda2.max)) lambda2.max <- 1e3*lambda2
  if(is.null(lambda1.max)) lambda1.max <- 1e3*lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2), length.out = length)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) alpha0 <- rep(0, q)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  mu <- mu0
  alpha <- alpha0
  alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
  theta <- theta0
  iter <- 1
  list.mu <- list()
  list.alpha <- list()
  list.theta <- list()
  list.mu[[1]] <- mu0
  list.alpha[[1]] <- alpha0
  list.theta[[1]] <- theta0
  list.res <- list()
  list.res[[1]] <- mu0 + alpha.mat + theta0
  for(i in 2:length){
    res <- rwls_l1_nuc.cov(y, x, var.type = var.type, lambda1 = exp(lambda1.grid.log[i]),
                      lambda2 = exp(lambda2.grid.log[i]), maxit = maxit,
                      upper = upper, lower = lower, mu0 = list.mu[[iter]], alpha0 = list.alpha[[iter]],
                      theta0 = list.theta[[iter]], thresh = thresh, trace.it = trace.it,
                      scale = scale, offset = offset)
    iter <- iter + 1
    list.res[[iter]] <- res
    list.mu[[iter]] <- res$mu
    list.alpha[[iter]] <- res$alpha
    list.theta[[iter]] <- res$theta
  }
  return(list(lambda1.grid = lambda1.grid.log, lambda2.grid = lambda2.grid.log,
              list.res = list.res, list.mu = list.mu, list.alpha = list.alpha,
              list.theta = list.theta))
}


#' mimi.multi
#' Compute solution of mimi with group effects along regularization path
#' @param y nxp observation matrix
#' @param groups factor of length n indicating group memberships
#' @param var.type vector of length p indicating types of y columns (gaussian, binomial, poisson)
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param lambda2 positive number, regularization parameter for l1 norm penalty
#' @param maxit integer, maximum number of iterations
#' @param mu0 real number initial offset (optional)
#' @param alpha0 matrix of size Nxp (N nb of groups) initial group effects (optional)
#' @param theta0 matrix of size nxp, initial interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean, whether convergence information should be printed
#' @param lambda1.max positive number, max value of regularization parameter for nuclear norm penalty
#' @param lambda2.max positive number, max value of regularization parameter for l1 norm penalty
#' @param length integer, size of grid for regularization path
#' @param upper real number, upper bound on parameters entries
#' @param lower real number, lower bound on parameters entries
#' @param scale boolean, whether or not cost functions should be scaled by column
#' @param offset boolean, whether or not an offset should be fitted, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 (log scale)}
#' \item{list.res}{List of same length as regularization grid with output of mimi for each lambda}
#' \item{list.mu}{List of same length as regularization grid with offset for each lambda}
#' \item{list.alpha}{List of same length as regularization grid with value of alpha for each lambda}
#' \item{list.theta}{List of same length as regularization grid with value of theta for each lambda}
#' @export
#'
#' @examples
#' y0 <- matrix(rnorm(6 * 10), nrow = 6)
#' y0[sample(1:50, size = 10)] <- NA
#' groups <- c(1,1,2,2,3,3)
#' groups <- as.factor(groups)
#' var.type <- rep("gaussian", 10)
#' res <- mimi.multi(y0, groups, var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi.multi <- function(y, groups, var.type, lambda1,
                       lambda2, maxit = 100, mu0 = NULL,
                       alpha0 = NULL, theta0 = NULL, thresh = 1e-4,
                       trace.it = F, lambda1.max = NULL,
                       lambda2.max = NULL, length = 20,
                       upper = 12, lower = -12, scale = F,
                       offset = F, max.rank = 5)
{
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y)
  w.groups <- aggregate(matrix(omega, nrow = n), list(groups), sum)[, 2:(p + 1)]
  w.groups[w.groups == 0] <- 1
  if(is.null(lambda2.max)) lambda2.max <- 1e3*lambda2
  ncenters <- aggregate(rep(1, n), list(groups), sum)[,2]
  N <- nlevels(groups)
  if(is.null(lambda1.max)) lambda1.max <- 1e3*lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2), length.out = length)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) {
    alpha0 <- rep(0, N * p)
    alpha.rep <- matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = p)
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
  list.mu <- list()
  list.alpha <- list()
  list.theta <- list()
  list.mu[[1]] <- mu0
  list.alpha[[1]] <- alpha0
  list.theta[[1]] <- theta0
  list.res <- list()
  list.res[[1]] <- mu0 + alpha.rep + theta0
  for(i in 2:length){
    res <- rwls_l1_nuc.multi(y, groups = groups, var.type = var.type,
                             lambda1 = exp(lambda1.grid.log[i]),
                             lambda2 = exp(lambda2.grid.log[i]),
                             maxit = maxit, upper = upper,
                             lower = lower, mu0 = list.mu[[iter]],
                             alpha0 = list.alpha[[iter]],
                             theta0 = list.theta[[iter]],
                             thresh = thresh, trace.it = trace.it,
                             scale = scale, offset = offset,
                             max.rank = 5)
    iter <- iter + 1
    list.res[[iter]] <- res
    list.mu[[iter]] <- res$mu
    list.alpha[[iter]] <- res$alpha
    list.theta[[iter]] <- res$theta
  }
  return(list(lambda1.grid = lambda1.grid.log, lambda2.grid = lambda2.grid.log,
              list.res = list.res, list.mu = list.mu, list.alpha = list.alpha,
              list.theta = list.theta))
}

