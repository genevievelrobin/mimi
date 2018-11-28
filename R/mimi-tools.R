

#' covmat
#'
#' @param R nxK1 matrix of row covariates
#' @param C nxK2 matrix of column covariates
#' @param n number of rows
#' @param p number ofcolumns
#'
#' @return the joint product of R and C, a (np)x(K1+K1) matrix in order row1col1,row2col1,...,rowncol1, row1col2, row2col2,...,rowncolp
#' @export
#'
#' @examples
#' R <- matrix(rnorm(10), 5)
#' C <- matrix(rnorm(9), 3)
#' covs <- covmat(R,C,5,3)
covmat <- function(R = NULL, C = NULL, n, p) {
  if (is.null(R)) {
    covs <- covmatC(C, n)
  } else if (is.null(C)) {
    covs <- covmatR(R, p)
  } else {
    R <- as.matrix(R)
    C <- as.matrix(C)
    dR <- dim(R)
    dC <- dim(C)
    K1 <- dR[2]
    K2 <- dC[2]
    covs <-
      cbind(do.call(rbind, replicate(nrow(C), R, simplify = FALSE)),
            C[rep(seq_len(nrow(C)), each = nrow(R)), ])
  }
  return(covs)
}
#' covmatR
#'
#' @param R nxK1 matrix of row covariates
#' @param p number ofcolumns
#'
#' @return repeats every row of R p times
#' @export
#'
#' @examples
#' R <- matrix(rnorm(10), 5)
#' cov <- covmatR(R,3)
covmatR <- function(R, p) {
  R <- as.matrix(R)
  dR <- dim(R)
  n <- dR[1]
  K1 <- dR[2]
  covs <- do.call(rbind, replicate(p, R, simplify = FALSE))
  return(covs)
}

#' covmatC
#'
#' @param C nxK2 matrix of column covariates
#' @param n number of rows
#'
#' @return repeats C n times
#' @export
#'
#' @examples
#' C <- matrix(rnorm(10), 5)
#' cov <- covmatC(C,3)
covmatC <- function(C, n) {
  C <- as.matrix(C)
  dC <- dim(C)
  p <- dC[1]
  K2 <- dC[2]
  covs <- C[rep(seq_len(p), each = n), ]
  return(covs)
}

