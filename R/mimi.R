#' mimi (Main effects and Interactions in Mixed and Incomplete data frames)
#' The method estimates main effects (group effects or effects of covariates)
#' and interactions in mixed data frames with missing values. The results can be used
#' for imputation or interpretation purposes.
#' @param y nxp matrix of observations
#' @param model either one of "groups", "covariates" or "low-rank", indicating which model should be fitted
#' @param x (np)xN matrix of covariates (optional)
#' @param groups factor of length n indicating groups (optional)
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binary or poisson)
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param maxit integer maximum number of iterations
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{alpha}{vector of main effects}
#' \item{theta}{interaction matrix}
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
#' idx_NA <- sample(1:(3 * n * p), size = round(0.01 * 3 * n * p))
#' y[idx_NA] <- NA
#' res <- mimi(y, model = "low-rank", var.type = var.type, lambda1 = 1, maxit=5)
mimi <-
  function(y,
           model = c("covariates", "low-rank"),
           x = NULL,
           groups = NULL,
           var.type = c("gaussian", "binary", "poisson"),
           lambda1,
           lambda2,
           maxit = 100,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    if (model == "covariates") {
      return(
        mimi.cov(
          y = y,
          x = x,
          var.type = var.type,
          lambda1 = lambda1,
          lambda2 = lambda2,
          maxit = maxit,
          alpha0 = alpha0,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
      )
    } else{
      return(
        mimi.lr(
          y = y,
          var.type = var.type,
          lambda1 = lambda1,
          maxit = maxit,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
      )

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
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{yimputed}{imputed data set}
#' \item{param}{estimated parameter matrix}
#' \item{alpha}{estimated vector of main effects}
#' \item{theta}{estimated interaction matrix}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
#' @examples
#' n = 6; p = 2
#' x <- matrix(rnorm(6*6*2), nrow = 6*6)
#' param <- matrix(x%*%c(2,0), nrow=6)+matrix(rnorm(6*6), nrow=6)
#' y1 <- matrix(rnorm(mean = c(param[, 1:2]), n * p), nrow = n)
#' y2 <- matrix(rbinom(n * p, prob = c(exp(param[,3:4])/(1+exp(param[,3:4]))), size = 1), nrow = n)
#' y3 <- matrix(rpois(n * p, lambda = c(exp(param[,5:6]))), nrow = n)
#' y <- cbind(y1, y2, y3)
#' var.type <- c(rep("gaussian", p), rep("binary", p), rep("poisson", p))
#' idx_NA <- sample(1:(3 * n * p), size = round(0.1 * 3 * n * p))
#' y[idx_NA] <- NA
#' res <- mimi.cov(y, x, var.type = var.type, lambda1 = 1, lambda2 = 2, thresh=0.1)
mimi.cov <-
  function(y,
           x,
           var.type = c("gaussian", "binary", "poisson"),
           lambda1,
           lambda2,
           maxit = 100,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    yy <- y
    n <- nrow(y)
    if (sum(var.type == "binary") > 0) {
      for (j in 1:sum(var.type == "binary")) {
        y <- data.frame(y)
        y[, which(var.type == "binary")[j]] <-
          as.factor(y[, which(var.type == "binary")[j]])
      }
    }
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n, p)-1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    q <- ncol(x)
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    alpha <- alpha0
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    theta <- theta0
    res <- irwls.cov(
      y,
      x,
      var.type = var.type,
      lambda1 = lambda1,
      lambda2 = lambda2,
      maxit = maxit,
      alpha0 = alpha,
      theta0 = theta,
      thresh = thresh,
      trace.it = trace.it
    )
    alpha <- res$alpha
    theta <- res$theta

    return(list(
      y.imputed = res$y.imputed,
      param = res$param,
      alpha = alpha,
      theta = theta
    ))
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
#' @param max.rank integer, maximum rank of interaction matrix
#' @return A list with the following elements
#' \item{yimputed}{the imputed data set}
#' \item{theta}{estimated low-rank matrix}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
#' @examples
#' y0 <- matrix(rnorm(6 * 10), nrow = 6)
#' y0[sample(1:50, size = 10)] <- NA
#' var.type <- rep("gaussian", 10)
#' res <- mimi.lr(y0, var.type, lambda1 = 0.1)
mimi.lr <-
  function(y,
           var.type = c("gaussian", "binary"),
           lambda1,
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    yy <- y
    if (sum(var.type == "binary") > 0) {
      for (j in 1:sum(var.type == "binary")) {
        y <- data.frame(y)
        y[, which(var.type == "binary")[j]] <-
          as.factor(y[, which(var.type == "binary")[j]])
      }
    }
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n,p)-1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    theta <- theta0
    res <-
      irwls.lr(
        y,
        var.type = var.type,
        lambda1 = lambda1,
        maxit = maxit,
        theta0 = theta,
        thresh = thresh,
        trace.it = trace.it,
        max.rank = max.rank
      )
    return(list(y.imputed = res$y.imputed, theta = res$theta))
  }


#' cv.mimi
#'
#' @param y [matrix, data.frame] incomplete and mixed data frame (nxp)
#' @param model either one of "groups", "covariates" or "low-rank", indicating which model should be fitted
#' @param var.type vector of length p indicating types of y columns (gaussian, binary, poisson)
#' @param x [matrix, data.frame] covariate matrix (npxq)
#' @param N [integer] number of cross-validation folds
#' @param thresh [positive number] convergence threshold, default is 1e-5
#' @param maxit [integer] maximum number of iterations, default is 100
#' @param max.rank [integer] maximum rank of interaction matrix, default is 2
#' @param trace.it [boolean] whether information about convergence should be printed
#' @param parallel [boolean] whether the N-fold cross-validation should be parallelized, default value is TRUE
#' @param len [integer] the size of the grid
#'
#' @return A list with the following elements
#' \item{lambda1}{regularization parameter estimated by cross-validation for nuclear norm penalty (interaction matrix)}
#' \item{lambda2}{regularization parameter estimated by cross-validation for l1 norm penalty (main effects)}
#' \item{errors}{a table containing the prediction errors for all pairs of parameters}

#' @export
#' @import data.table doParallel parallel foreach
cv.mimi <- function(y,
                    model=c("low-rank", "covariates"),
                    var.type,
                    x = NULL,
                    N = 10,
                    thresh = 1e-5,
                    maxit = 100,
                    max.rank = NULL,
                    trace.it = F,
                    parallel = T,
                    len = 20) {
  y <- as.matrix(y)
  Y2 <- y
  Y2[is.na(Y2)] <- 0
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  if(is.null(max.rank)) max.rank <- min(n,p)-1 else max.rank <- min(max.rank, min(n, p)-1)
  m <- sum(!is.na(y))
  na_func <- function(x, prob = 0.1) {
    x <- as.matrix(x)
    omega <- !is.na(x)
    obs.idx <- which(omega)
    yp <- x
    yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    return(yp)
  }
  lambda1.max <- max(svd(Y2)$d)
  lambda1.min <- 1e-3*lambda1.max
  grid.lambda1 <-
    exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
  if(model=="covariates"){
    lambda2.max <-
      max(y, na.rm=T)
    lambda2.min <- 1e-3 * lambda2.max
    grid.lambda2 <-
      exp(seq(log(lambda2.min), log(lambda2.max), length.out = len))
    grid <- as.matrix(data.table::CJ(grid.lambda1, grid.lambda2))
    grid <- grid[nrow(grid):1,]
    if (parallel) {
      nbco <- detectCores()
      cl <- makeCluster(nbco)
      registerDoParallel(cl)
      res.cv <-
        foreach(k = 1:N,
                .packages = c("mimi", "parallel", "glmnet")) %dopar% {
                  sapply(1:nrow(grid),
                         function(i) {
                           yy <- na_func(as.matrix(y), prob = 0.1)
                           if (trace.it)
                             print(paste("lambda", i))
                           res <-
                             mimi(as.matrix(yy),
                                  model="covariates",
                                  var.type=var.type,
                                  x=x,
                                  lambda1 = grid[i, 1],
                                  lambda2 = grid[i, 2])$y.imputed

                           return(sqrt(sum((res - y) ^ 2, na.rm = T)))
                         })
                }
    } else{
      ylist <-
        lapply(1:N, function(k)
          na_func(as.matrix(y), prob = 0.1))
      res.cv <-   lapply(1:N, function(k) {
        if (trace.it)
          print(paste("boot", k))
        sapply(1:nrow(grid),
               function(i) {
                 if (trace.it)
                   print(paste("lambda", i))
                 res <-
                   mimi(ylist[[k]], model="covariates", var.type=var.type, x=x, lambda1 = grid[i, 1], lambda2 = grid[i, 2])$y.imputed
                 return(sqrt(sum((res - y) ^ 2, na.rm = T)))
               })
      })
    }
    res.cv <- colMeans(do.call(rbind, res.cv))
    l <- which.min(res.cv)
    lambda <- grid[l, ]
    names(lambda) <- c("lambda1","lambda2")
    dat <-
      data.frame(errors = res.cv,
                 lambda1 = grid[, 1],
                 lambda2 = grid[, 2])
  } else{
    if (parallel) {
      nbco <- detectCores()
      cl <- makeCluster(nbco)
      registerDoParallel(cl)
      res.cv <-
        foreach(k = 1:N,
                .packages = c("mimi", "parallel", "glmnet")) %dopar% {
                  sapply(1:len,
                         function(i) {
                           yy <- na_func(as.matrix(y), prob = 0.1)
                           if (trace.it)
                             print(paste("lambda", i))
                           res <-
                             mimi(as.matrix(yy),
                                  model="low-rank", var.type=var.type,
                                  lambda1 = grid.lambda1[i])$y.imputed
                           return(sqrt(sum((res - y) ^ 2, na.rm = T)))
                         })
                }
    } else{
      ylist <-
        lapply(1:N, function(k)
          na_func(as.matrix(y), prob = 0.1))
      res.cv <-   lapply(1:N, function(k) {
        if (trace.it)
          print(paste("boot", k))
        sapply(1:len,
               function(i) {
                 if (trace.it)
                   print(paste("lambda", i))
                 res <-
                   mimi(ylist[[k]], model="low-rank", var.type=var.type, lambda1 = grid.lambda1[i])$y.imputed
                 return(sqrt(sum((res - y) ^ 2, na.rm = T)))
               })
      })
    }
    res.cv <- colMeans(do.call(rbind, res.cv))
    l <- which.min(res.cv)
    lambda <- grid.lambda1[l]
    dat <-
      data.frame(errors = res.cv,
                 lambda1 = grid.lambda1)

  }

  return(list(
    lambda = lambda,
    errors = dat
  ))
}
