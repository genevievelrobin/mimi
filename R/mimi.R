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
#' @param mu0 real number initial value of offset (optional)
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean indicating whether or not an offset should be fitted
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{mu}{offset}
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
           model = c("groups", "covariates", "low-rank"),
           x = NULL,
           groups = NULL,
           var.type = c("gaussian", "binary", "categorical", "poisson"),
           lambda1,
           lambda2,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-6,
           trace.it = F,
           offset = T,
           scale = T,
           max.rank = 5)
  {
    if (model == "groups") {
      return (
        mimi.multi(
          y = y,
          groups = groups,
          var.type = var.type,
          lambda1 = lambda1,
          lambda2 = lambda2,
          maxit = maxit,
          mu0 = mu0,
          alpha0 = alpha0,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          offset = offset,
          scale = scale,
          max.rank = max.rank
        )
      )
    } else if (model == "covariates") {
      return(
        mimi.cov(
          y = y,
          x = x,
          var.type = var.type,
          lambda1 = lambda1,
          lambda2 = lambda2,
          maxit = maxit,
          mu0 = mu0,
          alpha0 = alpha0,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          offset = offset,
          scale = scale,
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
          offset = offset,
          scale = scale,
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
#' @param mu0 real number initial value of offset (optional)
#' @param alpha0 vector of length N: initial value of regression parameter (optional)
#' @param theta0 matrix of size nxp: initial value of interactions (optional)
#' @param thresh positive number, convergence criterion
#' @param trace.it boolean indicating whether convergence information should be printed
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean indicating whether or not an offset should be fitted
#' @param max.rank integer, maximum rank of interaction matrix theta
#'
#' @return A list with the following elements
#' \item{yimputed}{imputed data set}
#' \item{param}{estimated parameter matrix}
#' \item{mu}{estimated offset}
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
           var.type = c("gaussian", "binary", "categorical", "poisson"),
           lambda1,
           lambda2,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-4,
           trace.it = F,
           offset = T,
           scale = F,
           max.rank = 5)
  {
    yy <- y
    nlevel = rep(1, ncol(y))
    n <- nrow(y)
    if (sum(var.type == "binary") > 0) {
      for (j in 1:sum(var.type == "binary")) {
        y <- data.frame(y)
        y[, which(var.type == "binary")[j]] <-
          as.factor(y[, which(var.type == "binary")[j]])
        y[, which(var.type == "binary")[j]] <-
          sapply(as.character(y[, which(var.type == "binary")[j]]), function(t)
            if (t %in% c(levels(y[, which(var.type == "binary")[j]])[1]))
              1
            else if (is.na(t))
              NA
            else
              0)
      }
    }
    if (sum(var.type == "categorical") > 0) {
      y <- data.frame(y)
      for (j in 1:sum(var.type == "categorical")) {
        y[, which(var.type == "categorical")[j]] <-
          as.factor(y[, which(var.type == "categorical")[j]])
      }
      tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
      nlevel[var.type == "categorical"] <-
        sapply(which(var.type == "categorical"), function(t)
          nlevels(y[, t]))
      colnames(tab) <-
        paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
              , ".", c(sapply(nlevel[var.type == "categorical"], function(k)
                1:k)), sep = "")
      y2 <- rep(0, n)
      x2 <- rep(0, ncol(x))
      vt2 <- NULL
      count <- 1
      for (j in 1:ncol(y)) {
        if (!(var.type[j] == "categorical")) {
          y2 <- cbind(y2, y[j])
          vt2 <- c(vt2, var.type[j])
          x2 <- rbind(x2, x[(1 + j - 1):(j + n - 1),])
        } else{
          y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) - 1)])
          vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
          for (i in 1:nlevels(y[, j])) {
            x2 <- rbind(x2, x[(1 + j - 1):(j + n - 1),])
          }
        }
      }
      x <- x2[2:nrow(x2),]
      y <- y2[, 2:ncol(y2)]
    } else
      vt2 <- var.type
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    max.rank <- min(max.rank, min(n, p) - 1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    q <- ncol(x)
    if (is.null(mu0))
      mu0 <- 0
    if (is.null(alpha0))
      alpha0 <- rep(0, q)
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    mu <- mu0
    alpha <- alpha0
    alpha.mat <-
      matrix(matrix(as.numeric(x), nrow = n * p) %*% alpha, nrow = n)
    theta <- theta0
    res <- irwls.cov(
      y,
      x,
      var.type = var.type,
      nlevel = nlevel,
      lambda1 = lambda1,
      lambda2 = lambda2,
      maxit = maxit,
      mu0 = mu,
      alpha0 = alpha,
      theta0 = theta,
      thresh = thresh,
      trace.it = trace.it,
      offset = offset,
      scale = scale,
      vt2 = vt2
    )
    mu <- res$mu
    alpha <- res$alpha
    theta <- res$theta

    return(list(
      y.imputed = res$y.imputed,
      param = res$param,
      mu = mu,
      alpha = alpha,
      theta = theta
    ))
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
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean, whether or not an offset should be fitted, default FALSE
#' @param max.rank integer, maximum rank of interaction matrix
#'
#' @return A list with the following elements
#' \item{yimputed}{imputed data set}
#' \item{param}{estimated parameter matrix}
#' \item{mu}{estimated offset}
#' \item{alpha}{estimated matrix of group effects}
#' \item{theta}{estimated interaction matrix}
#' @export
#' @importFrom FactoMineR tab.disjonctif.prop
#' @examples
#' param <- matrix(rnorm(6 * 10), nrow = 6)
#' groups <- c(1,1,2,2,3,3)
#' groups <- as.factor(groups)
#' param[groups==1, 1] <- 2
#' param[groups==2, 2] <- 2
#' y <- matrix(rnorm(6*10,c(param)), nrow=6)
#' var.type <- rep("gaussian", 10)
#' res <- mimi.multi(y, groups, var.type, lambda1 = 1, lambda2 = 2, thresh=0.1)
mimi.multi <-
  function(y,
           groups,
           var.type = c("gaussian", "binary", "categorical", "poisson"),
           lambda1,
           lambda2,
           maxit = 100,
           mu0 = NULL,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           offset = T,
           scale = T,
           max.rank = 5)
  {
    yy <- y
    n <- nrow(y)
    nlevel = rep(1, ncol(y))
    if (sum(var.type == "binary") > 0) {
      for (j in 1:sum(var.type == "binary")) {
        y <- data.frame(y)
        y[, which(var.type == "binary")[j]] <-
          as.factor(y[, which(var.type == "binary")[j]])
        y[, which(var.type == "binary")[j]] <-
          sapply(as.character(y[, which(var.type == "binary")[j]]), function(t)
            if (t %in% c(levels(y[, which(var.type == "binary")[j]])[1]))
              1
            else if (is.na(t))
              NA
            else
              0)
      }
    }
    if (sum(var.type == "categorical") > 0) {
      y <- data.frame(y)
      for (j in 1:sum(var.type == "categorical")) {
        y[, which(var.type == "categorical")[j]] <-
          as.factor(y[, which(var.type == "categorical")[j]])
      }
      tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
      tab[((tab < 1) * (tab > 0)) == 1] <- NA
      nlevel[var.type == "categorical"] <-
        sapply(which(var.type == "categorical"), function(t)
          nlevels(y[, t]))
      colnames(tab) <-
        paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
              , ".", c(sapply(nlevel[var.type == "categorical"], function(k)
                1:k)), sep = "")
      vt2 <- NULL
      y2 <- rep(0, n)
      count <- 1
      for (j in 1:ncol(y)) {
        if (!(var.type[j] == "categorical")) {
          y2 <- cbind(y2, y[j])
          vt2 <- c(vt2, var.type[j])
        } else{
          y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) - 1)])
          vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
          count <- count + nlevels(y[, j])
        }
      }
      y <- y2[, 2:ncol(y2)]
    } else
      vt2 <- var.type
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    max.rank <- min(max.rank, min(n, p) - 1)
    y <- as.matrix(y)
    y <- matrix(as.numeric(y), nrow = n)
    omega <- !is.na(y)
    ncenters <- aggregate(rep(1, n), list(groups), sum)[, 2]
    N <- nlevels(groups)
    if (is.null(mu0))
      mu0 <- 0
    if (is.null(alpha0)) {
      alpha.rep <- matrix(rep(0, n * p), nrow = n)
      alpha0 <- matrix(rep(0, N * p), nrow = N)
    } else{
      alpha.rep <-
        matrix(rep(as.matrix(alpha0), rep(ncenters, p)), nrow = n)
      alpha <- matrix(as.matrix(alpha0), nrow = N)
    }
    if (is.null(theta0))
      theta0 <- matrix(rep(0, n * p), nrow = n)
    mu <- mu0
    alpha <- alpha0
    alpha.rep <-
      matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
    theta <- theta0
    res <-
      irwls.multi(
        y,
        groups = groups,
        var.type = var.type,
        nlevel = nlevel,
        lambda1 = lambda1,
        lambda2 = lambda2,
        maxit = maxit,
        mu0 = mu,
        alpha0 = alpha,
        theta0 = theta,
        thresh = thresh,
        trace.it = trace.it,
        offset = offset,
        max.rank = max.rank,
        vt2 = vt2,
        scale = scale
      )
    mu <- res$mu
    alpha <- res$alpha
    theta <- res$theta
    return(list(
      y.imputed = res$y.imputed,
      param = res$param,
      mu = mu,
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
#' @param scale boolean indicating whether or not the column loss functions should be scaled
#' @param offset boolean, whether or not an offset should be fitted, default FALSE
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
           var.type = c("gaussian", "binary", "categorical", "poisson"),
           lambda1,
           maxit = 100,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           offset = T,
           scale = T,
           max.rank = 5)
  {
    yy <- y
    nlevel = rep(1, ncol(y))
    if (sum(var.type == "binary") > 0) {
      for (j in 1:sum(var.type == "binary")) {
        y <- data.frame(y)
        y[, which(var.type == "binary")[j]] <-
          as.factor(y[, which(var.type == "binary")[j]])
        y[, which(var.type == "binary")[j]] <-
          sapply(as.character(y[, which(var.type == "binary")[j]]), function(t)
            if (t %in% c(levels(y[, which(var.type == "binary")[j]])[1]))
              1
            else if (is.na(t))
              NA
            else
              0)
      }
    }
    if (sum(var.type == "categorical") > 0) {
      y <- data.frame(y)
      for (j in 1:sum(var.type == "categorical")) {
        y[, which(var.type == "categorical")[j]] <-
          as.factor(y[, which(var.type == "categorical")[j]])
      }
      tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
      tab[((tab > 0) * (tab < 1) == 1)] <- NA
      nlevel[var.type == "categorical"] <-
        sapply(which(var.type == "categorical"), function(t)
          nlevels(y[, t]))
      colnames(tab) <-
        paste(rep(colnames(y[, var.type == "categorical"]), nlevel[var.type == "categorical"])
              , ".", c(sapply(nlevel[var.type == "categorical"], function(k)
                1:k)), sep = "")
      y2 <- rep(0, n)
      vt2 <- NULL
      count <- 1
      for (j in 1:ncol(y)) {
        if (!(var.type[j] == "categorical")) {
          y2 <- cbind(y2, y[j])
          vt2 <- c(vt2, var.type[j])
        } else{
          y2 <- cbind(y2, tab[, count:(count + nlevels(y[, j]) - 1)])
          vt2 <- c(vt2, rep("categorical", nlevels(y[, j])))
        }
      }
      y <- y2[, 2:ncol(y2)]
    } else
      vt2 <- var.type
    y <- as.matrix(y)
    d <- dim(y)
    n <- d[1]
    p <- d[2]
    max.rank <- min(max.rank, min(n, p) - 1)
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
        nlevel = nlevel,
        lambda1 = lambda1,
        maxit = maxit,
        theta0 = theta,
        thresh = thresh,
        trace.it = trace.it,
        offset = offset,
        max.rank = max.rank,
        vt2 = vt2,
        scale = scale
      )
    theta <- res$theta

    return(list(y.imputed = res$y.imputed, theta = theta))
  }
