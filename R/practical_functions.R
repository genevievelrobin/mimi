#' covmat
#'
#' @param XR (nxqr) matrix of row covariates
#' @param XC (pxqc) matrix of column covariates
#'
#' @return X, a ((np)x(qr+qc)) matrix combining XR and XC
#' @export
#'
#' @examples
#' XR <- matrix(rnorm(10*3), nrow = 10)
#' XC <- matrix(rnorm(15*2), nrow = 15)
#' X <- covmat(XR, XC)
covmat <- function(XR = NULL, XC = NULL){
  XR <- as.matrix(XR)
  XC <- as.matrix(XC)
  dr <- dim(XR)
  n <- dr[1]
  qr <- dr[2]
  dc <- dim(XC)
  p <- dc[1]
  qc <- dc[2]
  XRrep <- matrix(rep(XR, p), nrow = p*n, byrow = T)
  XCrep <- matrix(rep(XC, each = n), nrow = p*n)
  return(cbind(XRrep, XCrep))
}
