#' main function: low-rank models to analyze and impute mixed and incomplete data frames
#' with numeric, binary and discrete variables, and missing values
#' @param y nxp matrix of observations
#' @param model either one of "groups", "covariates" or "low-rank", indicating which model should be fitted
#' @param x (np)xN matrix of covariates (optional)
#' @param groups factor of length n indicating groups (optional)
#' @param var.type vector of length p indicating the data types of the columns of y (gaussian, binomial or poisson)
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param algo type of algorithm to use, either one of "bcgd" (small dimensions, gaussian and binomial variables) or "mcgd" (large dimensions, poisson variables)
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
#' var.type <- c(rep("gaussian", p), rep("binomial", p), rep("poisson", p))
#' idx_NA <- sample(1:(3 * n * p), size = round(0.01 * 3 * n * p))
#' y[idx_NA] <- NA
#' res <- mimi(y, model = "low-rank", var.type = var.type, lambda1 = 1, maxit=5)
mimi <-
  function(y,
           model = c("low-rank", "multilevel", "covariates"),
           x = NULL,
           groups = NULL,
           var.type = c("gaussian", "binomial", "poisson"),
           lambda1,
           lambda2,
           algo = c("mcgd", "bcgd"),
           maxit = 100,
           alpha0 = NULL,
           theta0 = NULL,
           thresh = 1e-5,
           trace.it = F,
           max.rank = NULL)
  {
    model <- match.arg(model,c("low-rank","multilevel", "covariates"),several.ok=T)[1]
    if (model == "covariates") {
      return(
        mimi.cov(
          y = y,
          x = x,
          var.type = var.type,
          lambda1 = lambda1,
          lambda2 = lambda2,
          algo=algo,
          maxit = maxit,
          alpha0 = alpha0,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
      )
    } else if(model=="low-rank"){
      return(
        mimi.lr(
          y = y,
          var.type = var.type,
          lambda1 = lambda1,
          algo=algo,
          maxit = maxit,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
      )

    } else if(model == "multilevel"){
      return(
        mimi.ml(
          y = y,
          groups=groups,
          var.type = var.type,
          lambda1 = lambda1,
          lambda2 = lambda2,
          algo=algo,
          maxit = maxit,
          theta0 = theta0,
          thresh = thresh,
          trace.it = trace.it,
          max.rank = max.rank
        )
      )

    }
  }


