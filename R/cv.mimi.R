#' selection of the regularization
#' parameters (lambda1 and lambda2) of the mimi function by cross-validation
#'
#' @param y [matrix, data.frame] incomplete and mixed data frame (nxp)
#' @param model either one of "groups", "covariates" or "low-rank", indicating which model should be fitted
#' @param var.type vector of length p indicating types of y columns (gaussian, binomial, poisson)
#' @param x [matrix, data.frame] covariate matrix (npxq)
#' @param groups factor of length n indicating groups (optional)
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
                    model = c("low-rank", "covariates"),
                    var.type,
                    x = NULL,
                    groups = NULL,
                    N = 5,
                    thresh = 1e-5,
                    maxit = 100,
                    max.rank = NULL,
                    trace.it = F,
                    parallel = T,
                    len = 15) {
  y <- as.matrix(y)
  Y2 <- y
  Y2[is.na(Y2)] <- 0
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  if (is.null(max.rank))
    max.rank <- min(n, p) - 1
  else
    max.rank <- min(max.rank, min(n, p) - 1)
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
  lambda1.min <- 1e-3 * lambda1.max
  grid.lambda1 <-
    exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
  if ((model == "covariates") || (model == "multilevel")) {
    if(model=="covariates"){
      gd_main <- grad_alpha(Y, cov, alpha, theta, var.type)
      lambda2.max <- max(abs(gd_main))
    } else{
      lambda2.max <- max(aggregate(y, by=list(groups), sum, na.rm=T)[,2:(p+1)], na.rm = T)
    }
    lambda2.min <- max(1e-4, 1e-3 * lambda2.max)
    grid.lambda2 <-
      exp(seq(log(lambda2.min), log(lambda2.max), length.out = len))
    grid <- as.matrix(data.table::CJ(grid.lambda1, grid.lambda2))
    grid <- grid[nrow(grid):1, ]
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
                             mimi(
                               as.matrix(yy),
                               model = model,
                               var.type = var.type,
                               x = x,
                               groups = groups,
                               lambda1 = grid[i, 1],
                               lambda2 = grid[i, 2]
                             )$y.imputed

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
                   mimi(
                     ylist[[k]],
                     model = model,
                     var.type = var.type,
                     x = x,
                     groups=groups,
                     lambda1 = grid[i, 1],
                     lambda2 = grid[i, 2]
                   )$y.imputed
                 return(sqrt(sum((res - y) ^ 2, na.rm = T)))
               })
      })
    }
    res.cv <- colMeans(do.call(rbind, res.cv))
    l <- which.min(res.cv)
    lambda <- grid[l,]
    names(lambda) <- c("lambda1", "lambda2")
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
                             mimi(
                               as.matrix(yy),
                               model = "low-rank",
                               var.type = var.type,
                               lambda1 = grid.lambda1[i]
                             )$y.imputed
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
                   mimi(
                     ylist[[k]],
                     model = "low-rank",
                     var.type = var.type,
                     lambda1 = grid.lambda1[i]
                   )$y.imputed
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

  return(list(lambda = lambda,
              errors = dat))
}
