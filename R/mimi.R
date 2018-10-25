#' mimi
#' Compute solution of mimi with covariates/groups along regularization path
#' @param y nxp matrix of observations
#' @param model either one of "groups" or "covariates", indicating which model should be fitted
#' @param x (np)xN matrix of covariates (optional)
#' @param groups factor of length n indicating groups (optional)
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binary or poisson)
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
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean indicating whether or not an offset should be fitted
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 used for warm start (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 used for warm start (log scale)}
#' \item{list.res}{list of results of rwls_l1_nuc (same length as grid of lambda)}
#' \item{list.mu}{list of offsets (same length as grid of lambda)}
#' \item{list.alpha}{list of regression parameters (same length as grid of lambda)}
#' \item{list.theta}{list of interaction matrices (same length as grid of lambda)}
#' @export
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
#' res <- mimi(y, model = "covariates", x = x, var.type = var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi <- function(y, model = c("groups", "covariates", "low-rank"), x = NULL, groups = NULL,
                 var.type = c("gaussian", "binary", "categorical", "poisson"),
                 lambda1, lambda2, maxit = 100, mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                 thresh = 1e-6, trace.it = F, lambda1.max = NULL,
                 lambda2.max = NULL, length = 20,  offset = T, scale = T,
                 max.rank = 20, wmax = NULL)
{
  if(model == "groups"){
    return (mimi.multi(y = y, groups = groups, var.type = var.type, lambda1 = lambda1,
                       lambda2 = lambda2, maxit = maxit, mu0 = mu0, alpha0 = alpha0,
                       theta0 = theta0, thresh = thresh, trace.it = trace.it,
                       lambda1.max = lambda1.max, lambda2.max = lambda2.max,
                       length = length, offset = offset, scale = scale,
                       max.rank = max.rank, wmax = wmax))
  } else if(model == "covariates"){
    return(mimi.cov(y = y, x = x, var.type = var.type, lambda1 = lambda1,
                    lambda2 = lambda2, maxit = maxit, mu0 = mu0, alpha0 = alpha0,
                    theta0 = theta0, thresh = thresh, trace.it = trace.it,
                    lambda1.max = lambda1.max, lambda2.max = lambda2.max,
                    length = length, offset = offset, scale = scale,
                    max.rank = max.rank, wmax = wmax))
  } else{
    return(mimi.lr(y = y, var.type = var.type, lambda1 = lambda1, maxit = maxit,
                  theta0 = theta0, thresh = thresh, trace.it = trace.it,
                  lambda1.max = lambda1.max, length = length, offset = offset,
                  scale = scale, max.rank = max.rank, wmax = wmax))

  }
}



#' mimi.cov
#' Compute solution of mimi with covariates effects along regularization path
#' @param y nxp matrix of observations
#' @param x (np)xN matrix of covariates
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binary or poisson)
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
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean indicating whether or not an offset should be fitted
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 used for warm start (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 used for warm start (log scale)}
#' \item{list.res}{list of results of rwls_l1_nuc (same length as grid of lambda)}
#' \item{list.mu}{list of offsets (same length as grid of lambda)}
#' \item{list.alpha}{list of regression parameters (same length as grid of lambda)}
#' \item{list.theta}{list of interaction matrices (same length as grid of lambda)}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
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
#' res <- mimi.cov(y, x, var.type = var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi.cov <- function(y, x, var.type = c("gaussian", "binary", "categorical", "poisson"),
                     lambda1, lambda2, maxit = 1e3, mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                     thresh = 1e-4, trace.it = F, lambda1.max = NULL, lambda2.max = NULL,
                     length = 50, offset = T, scale = F, max.rank = 20, wmax = NULL)
{
  yy <- y
  nlevel = rep(1, ncol(y))
  n <- nrow(y)
  if(sum(var.type == "binary")>0){
    for(j in 1:sum(var.type == "binary")){
      y <- data.frame(y)
      y[, which(var.type == "binary")[j]] <- as.factor(y[, which(var.type == "binary")[j]])
      y[, which(var.type == "binary")[j]] <- sapply(as.character(y[, which(var.type == "binary")[j]]), function(t) if (t%in%c(levels(y[, which(var.type == "binary")[j]])[1])) 1 else if(is.na(t)) NA else 0)
    }
  }
  if(sum(var.type == "categorical")>0){
    y <- data.frame(y)
    for(j in 1:sum(var.type == "categorical")){
      y[, which(var.type == "categorical")[j]] <- as.factor(y[, which(var.type == "categorical")[j]])
    }
    tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
    nlevel[var.type == "categorical"] <- sapply(which(var.type == "categorical"), function(t) nlevels(y[, t]))
    colnames(tab) <- paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
                           ,".", c(sapply(nlevel[var.type == "categorical"], function(k) 1:k)), sep = "")
    y2 <- rep(0, n)
    x2 <- rep(0, ncol(x))
    vt2 <- NULL
    count <- 1
    for(j in 1:ncol(y)){
      if(!(var.type[j]=="categorical")){
        y2 <- cbind(y2, y[j])
        vt2 <- c(vt2, var.type[j])
        x2 <- rbind(x2, x[(1+j-1):(j+n-1), ])
      } else{
        y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) -1)])
        vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
        for(i in 1:nlevels(y[, j])){
          x2 <- rbind(x2, x[(1+j-1):(j+n-1), ])
        }
      }
    }
    x <- x2[2:nrow(x2), ]
    y <- y2[, 2:ncol(y2)]
  } else vt2 <- var.type
  y <- as.matrix(y)
  if(is.null(wmax)) wmax = 2*max(y, na.rm = T)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  max.rank <- min(max.rank, min(n,p)-1)
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  omega <- !is.na(y)
  q <- ncol(x)
  if(is.null(lambda2.max)) lambda2.max <- 1e4*lambda2
  if(is.null(lambda1.max)) lambda1.max <- 1e4*lambda1
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
    res <- rwls_l1_nuc.cov(y, x, var.type = var.type, nlevel = nlevel,
                           lambda1 = exp(lambda1.grid.log[i]),
                           lambda2 = exp(lambda2.grid.log[i]), maxit = maxit,
                           mu0 = list.mu[[iter]],
                           alpha0 = list.alpha[[iter]],
                           theta0 = list.theta[[iter]], thresh = thresh,
                           trace.it = trace.it, offset = offset, scale = scale, vt2 = vt2,
                           wmax = wmax)
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
#' @param var.type vector of length p indicating types of y columns (gaussian, binary, poisson)
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
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean, whether or not an offset should be fitted, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix

#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 (log scale)}
#' \item{lambda2.grid}{the grid of lambda2 (log scale)}
#' \item{list.res}{List of same length as regularization grid with output of mimi for each lambda}
#' \item{list.mu}{List of same length as regularization grid with offset for each lambda}
#' \item{list.alpha}{List of same length as regularization grid with value of alpha for each lambda}
#' \item{list.theta}{List of same length as regularization grid with value of theta for each lambda}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
#' @examples
#' y0 <- matrix(rnorm(6 * 10), nrow = 6)
#' y0[sample(1:50, size = 10)] <- NA
#' groups <- c(1,1,2,2,3,3)
#' groups <- as.factor(groups)
#' var.type <- rep("gaussian", 10)
#' res <- mimi.multi(y0, groups, var.type, lambda1 = 0.1, lambda2 = 0.2)
mimi.multi <- function(y, groups, var.type = c("gaussian", "binary", "categorical", "poisson"),
                       lambda1, lambda2, maxit = 100, mu0 = NULL, alpha0 = NULL, theta0 = NULL,
                       thresh = 1e-5, trace.it = F, lambda1.max = NULL, lambda2.max = NULL,
                       length = 20, offset = T, scale = T, max.rank = 20,
                       wmax = NULL)
{
  yy <- y
  n <- nrow(y)
  nlevel = rep(1, ncol(y))
  if(sum(var.type == "binary")>0){
    for(j in 1:sum(var.type == "binary")){
      y <- data.frame(y)
      y[, which(var.type == "binary")[j]] <- as.factor(y[, which(var.type == "binary")[j]])
      y[, which(var.type == "binary")[j]] <- sapply(as.character(y[, which(var.type == "binary")[j]]), function(t) if (t%in%c(levels(y[, which(var.type == "binary")[j]])[1])) 1 else if(is.na(t)) NA else 0)
    }
  }
  if(sum(var.type == "categorical")>0){
    y <- data.frame(y)
    for(j in 1:sum(var.type == "categorical")){
      y[, which(var.type == "categorical")[j]] <- as.factor(y[, which(var.type == "categorical")[j]])
    }
    tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
    tab[((tab<1)*(tab>0))==1] <- NA
    nlevel[var.type == "categorical"] <- sapply(which(var.type == "categorical"), function(t) nlevels(y[, t]))
    colnames(tab) <- paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
                           ,".", c(sapply(nlevel[var.type == "categorical"], function(k) 1:k)), sep = "")
    vt2 <- NULL
    y2 <- rep(0, n)
    count <- 1
    for(j in 1:ncol(y)){
      if(!(var.type[j]=="categorical")){
        y2 <- cbind(y2, y[j])
        vt2 <- c(vt2, var.type[j])
      } else{
        y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) -1)])
        vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
        count <- count + nlevels(y[, j])
      }
    }
    y <-y2[, 2:ncol(y2)]
  } else vt2 <- var.type
  #colnames(y) <- vt2
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  max.rank <- min(max.rank, min(n,p)-1)
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  if(is.null(wmax)) wmax = 2*max(y, na.rm = T)
  omega <- !is.na(y)
  if(is.null(lambda2.max)) lambda2.max <- 1e3*lambda2
  ncenters <- aggregate(rep(1, n), list(groups), sum)[,2]
  N <- nlevels(groups)
  if(is.null(lambda1.max)) lambda1.max <- 1e3*lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2), length.out = length)
  if(is.null(mu0)) mu0 <-0
  if(is.null(alpha0)) {
    alpha.rep <- matrix(rep(0, n*p), nrow = n)
    alpha0 <- matrix(rep(0, N * p), nrow = N)
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
    res <- rwls_l1_nuc.multi(y, groups = groups, var.type = var.type, nlevel = nlevel,
                             lambda1 = exp(lambda1.grid.log[i]), lambda2 = exp(lambda2.grid.log[i]),
                             maxit = maxit, mu0 = list.mu[[iter]],
                             alpha0 = list.alpha[[iter]], theta0 = list.theta[[iter]],
                             thresh = thresh, trace.it = trace.it, offset = offset, max.rank = max.rank,
                             vt2 = vt2, scale = scale, wmax = wmax)
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

#' mimi.lr
#' Compute solution of mimi for low-rank model along regularization path
#' @param y nxp observation matrix
#' @param var.type vector of length p indicating types of y columns (gaussian, binary, poisson)
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param maxit integer, maximum number of iterations
#' @param theta0 matrix of size nxp, initial interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean, whether convergence information should be printed
#' @param lambda1.max positive number, max value of regularization parameter for nuclear norm penalty
#' @param length integer, size of grid for regularization path
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean, whether or not an offset should be fitted, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix

#' @return A list with the following elements
#' \item{lambda1.grid}{the grid of lambda1 (log scale)}
#' \item{list.res}{List of same length as regularization grid with output of mimi for each lambda}
#' \item{list.theta}{List of same length as regularization grid with value of theta for each lambda}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
#' @examples
#' y0 <- matrix(rnorm(6 * 10), nrow = 6)
#' y0[sample(1:50, size = 10)] <- NA
#' var.type <- rep("gaussian", 10)
#' res <- mimi.lr(y0, var.type, lambda1 = 0.1)
mimi.lr <- function(y, var.type = c("gaussian", "binary", "categorical", "poisson"),
                       lambda1, maxit = 100, theta0 = NULL, thresh = 1e-5,
                       trace.it = F, lambda1.max = NULL, length = 20,
                       offset = T, scale = T, max.rank = 20, wmax = NULL)
{
  yy <- y
  nlevel = rep(1, ncol(y))
  if(sum(var.type == "binary")>0){
    for(j in 1:sum(var.type == "binary")){
      y <- data.frame(y)
      y[, which(var.type == "binary")[j]] <- as.factor(y[, which(var.type == "binary")[j]])
      y[, which(var.type == "binary")[j]] <- sapply(as.character(y[, which(var.type == "binary")[j]]), function(t) if (t%in%c(levels(y[, which(var.type == "binary")[j]])[1])) 1 else if(is.na(t)) NA else 0)
    }
  }
  if(sum(var.type == "categorical")>0){
    y <- data.frame(y)
    for(j in 1:sum(var.type == "categorical")){
      y[, which(var.type == "categorical")[j]] <- as.factor(y[, which(var.type == "categorical")[j]])
    }
    tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
    tab[((tab>0)*(tab<1) == 1)] <- NA
    nlevel[var.type == "categorical"] <- sapply(which(var.type == "categorical"), function(t) nlevels(y[, t]))
    colnames(tab) <- paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
                           ,".", c(sapply(nlevel[var.type == "categorical"], function(k) 1:k)), sep = "")
    y2 <- rep(0, n)
    x2 <- rep(0, ncol(x))
    vt2 <- NULL
    count <- 1
    for(j in 1:ncol(y)){
      if(!(var.type[j]=="categorical")){
        y2 <- cbind(y2, y[j])
        vt2 <- c(vt2, var.type[j])
      } else{
        y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) -1)])
        vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
      }
    }
    y <- y2[, 2:ncol(y2)]
  } else vt2 <- var.type
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  max.rank <- min(max.rank, min(n,p)-1)
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  if(is.null(wmax)) wmax = 2*max(y, na.rm = T)
  omega <- !is.na(y)
  if(is.null(lambda1.max)) lambda1.max <- 1e3*lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  if(is.null(theta0)) theta0 <- matrix(rep(0, n * p), nrow = n)
  theta <- theta0
  iter <- 1
  list.theta <- list()
  list.theta[[1]] <- theta0
  list.res <- list()
  list.res[[1]] <- list(y = y, theta = theta0, objective = 0, y.imputed = apply(y, c(1,2), function(xx){if(is.na(xx)) 0 else xx}))
  for(i in 2:length){
    res <- rwls_l1_nuc.lr(y, var.type = var.type, nlevel = nlevel, lambda1 = exp(lambda1.grid.log[i]),
                          maxit = maxit, theta0 = list.theta[[iter]], thresh = thresh,
                          trace.it = trace.it, offset = offset, max.rank = max.rank,
                          vt2 = vt2, scale = scale, wmax = wmax)
    iter <- iter + 1
    list.res[[iter]] <- res
    list.theta[[iter]] <- res$theta
  }
  return(list(lambda1.grid = lambda1.grid.log, list.res = list.res, list.theta = list.theta))
}

