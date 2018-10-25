#' wght
#'
#' @param y real number: observation
#' @param param real number: current value of parameter
#' @param var.type type of variable y (gaussian, binary, poisson)
#' @param nlevel vector indicating, for every column of y, the number of parameters of the distribution (1 for numerics, K for categorical with K categories)
#' @return weight of the quadratic approximation
#' @export
#' @examples
#' y <- rnorm(1)
#' param <- 0.5
#' var.type = "gaussian"
#' wght <- wght(y, param, var.type, nlevel = 1)
wght <- function(y, param, var.type, nlevel = NULL, wmax){
  if(var.type == "gaussian"){
    ytilde <- y
    vtilde2 <- 1
  } else if(var.type == "binary"){
    vtilde2 <- 1/4
    ytilde <- y/vtilde2 - (exp(param) / (1 + exp(param)))/vtilde2 + param
  } else if (var.type == "poisson"){
    vtilde2 <- exp(param)
    ytilde <- (y - exp(param))/ vtilde2  + param
  } else if(var.type == "categorical"){
    vtilde2 <- rep(0.25, nlevel)
    ytilde <- y/vtilde2 - (exp(param)/sum(exp(param)))/vtilde2 + param
  }
  else {
    print(var.type)
    stop("Incorrect type of variable. Should be 'gaussian', 'binary', 'poisson' or 'categorical'.")
  }
  return(list(ytilde = ytilde, vtilde2 = vtilde2))
}

#' quad_approx
#'
#' @param y observation matrix (nxP) - extended with dummies for cat variables
#' @param param matrix of parameters (nxP) - extended for cat variables
#' @param var.type type of the variables in y (length p<=P) - NOT extended for cat variables
#' @param nlevel vector indicating, for every column of y, the number of parameters of the distribution (1 for numerics, K for categorical with K categories) - NOT extended for cat variables
#' @return matrix of weights for the quadratic approximation (nxP)
#' @export
#' @examples
#' y <- matrix(rnorm(6*10), nrow = 6)
#' param <- matrix(rnorm(6*10), nrow = 6)
#' var.type <- rep("gaussian", 10)
#' nlevel <- rep(1, 10)
#' wghts <- quad_approx(y, param, var.type, nlevel)
quad_approx <- function(y, param, var.type, nlevel, wmax){
  n <- nrow(y)
  p <- length(var.type)
  yt <- rep(0, n)
  vt <- rep(0, n)
  count <- 1
  for (j in 1:p){
    w <- lapply(1:n, function(i) wght(y = y[i, count:(count+nlevel[j]-1)], param = param[i, count:(count+nlevel[j]-1)],
                                      var.type = var.type[j], nlevel = nlevel[j], wmax))
    count <- count + nlevel[j]
    if(nlevel[j] == 1){
      ytilde <- sapply(1:n, function(i) w[[i]]$ytilde)
      vtilde2 <- sapply(1:n, function(i) w[[i]]$vtilde2)
    } else{
      ytilde <- t(sapply(1:n, function(i) w[[i]]$ytilde))
      vtilde2 <- t(sapply(1:n, function(i) w[[i]]$vtilde2))
    }
    yt <- cbind(yt, as.matrix(ytilde))
    vt <- cbind(vt, as.matrix(vtilde2))
  }
  return(list(ytilde = yt[, 2:ncol(yt)], vtilde2 = vt[, 2:ncol(vt)]))
}

#' Frob
#' compute stopping criterion of weighted softImpute
#' @param Uold nxq matrix of left singular vectors
#' @param Dsqold vector singular values (q)
#' @param Vold pxq matrix of right singular vectors
#' @param U nxq matrix of left singular vectors
#' @param Dsq vector of singular values (q)
#' @param V pxq matrix of right singular vectors
#' @return value of stopping criterion
#' @export
#' @examples
#' U <- matrix(rnorm(5*2), nrow = 5)
#' V <- matrix(rnorm(3*2), nrow = 3)
#' Dsq <- abs(rnorm(2))
#' Uold <- matrix(rnorm(5*2), nrow = 5)
#' Vold <- matrix(rnorm(3*2), nrow = 3)
#' Dsqold <- abs(rnorm(2))
#' crit <- Frob(Uold,Dsqold,Vold,U,Dsq,V)
Frob <- function(Uold,Dsqold,Vold,U,Dsq,V){
  denom=sum(Dsqold^2)
  utu=Dsq* (t(U)%*%Uold)
  vtv=Dsqold* (t(Vold)%*%V)
  uvprod= sum(diag(utu%*%vtv))
  num=denom+sum(Dsq^2) -2*uvprod
  return(num/max(denom,1e-9))
}

#' log_factorial
#'
#' @param x matrix of integers
#' @return the log of the factorial of every entry in x
#' @export
#' @examples
#' x <- matrix(seq(1, 10), nrow = 5)
#' res <- log_factorial(x)
log_factorial <- function(x){
  m <- apply(x, c(1,2), function(t) if(is.na(t)){
    NA
  } else if(t==0){
    1
  } else  sum(log(1:t)))
  return(m)
}

#' wlra
#'
#' @param x nxp matrix to be approximates
#' @param w nxp matrix of weights (optional)
#' @param lambda regularization parameter
#' @param x0 nxp matrix of initial values
#' @param thresh convergence criterion
#' @param maxit maximum number of iterations
#' @param rank.max maximum number of sv to compute
#' @return a list containing the following elements
#' \item{d}{the vector of singular values}
#' \item{u}{the left singular vectors}
#' \item{v}{the right singular vectors}
#' \item{convergence}{boolean indicating if algorithm converged before maxit iterations}
#' \item{iter}{number of iterations performed}
#' @export
#' @examples
#' x <- matrix(rnorm(10*6), nrow = 10)
#' x[sample(60, 10)] <- NA
#' res <- wlra(x)
wlra <- function(x, w = NULL, lambda = 0, x0 = NULL, thresh = 1e-5, maxit = 1e3, rank.max = NULL)
  {
  d <- dim(x)
  n <- d[1]
  p <- d[2]
  if(is.null(w)) w = matrix(rep(1, n*p), nrow = n)
  if(is.null(rank.max)) rank.max <- min(n, p) - 1 else rank.max <- min(rank.max, min(n, p) - 1)
  if(is.null(x0)) x0 <- matrix(rep(0, n*p), nrow = n)
  xnas <- is.na(x)
  omega <- 1*(!xnas)
  nz=n*p-sum(xnas)
  xfill <- x
  xfill[xnas] <- 0
  xfill <- omega*w*xfill + (1 - omega*w)*x0
  xhat <- x0
  x0 <- xfill
  iter <- 0
  error <- 100
  svd.xfill=svd(xfill)
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    svd.old=svd(xhat)
    xhat.old <- xhat
    d=svd.xfill$d
    d=pmax(d-lambda,0)
    J <- min(rank.max, length(d))
    xhat <- svd.xfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.xfill$v[,seq(J)]))
    xfill <- omega*w*x0 + (1 - omega*w)*xhat
    svd.xfill=svd(xfill)
    denom <-  sum(w*(x-xhat.old)^2, na.rm = T) + lambda*sum(svd.old$d[seq(J)])
    error <- abs(sum(w*(x-xhat)^2, na.rm = T) + lambda*sum(d[seq(J)]) - denom)/denom
    #error=Frob(svd.old$u[, seq(J)],d[seq(J)],svd.old$v[, seq(J)],
    #           svd.xfill$u[, seq(J)],pmax(svd.xfill$d-lambda,0)[seq(J)],svd.xfill$v[, seq(J)])
  }
  if(error < thresh) cvg = T else cvg = F
  d <- svd.xfill$d[seq(J)]
  d=pmax(svd.xfill$d[seq(J)]-lambda,0)
  J=min(sum(d>0)+1,J)
  svd.xfill=list(u=svd.xfill$u[, seq(J)], d=d[seq(J)], v=svd.xfill$v[,seq(J)])
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  return(list(d = svd.xfill$d, u = svd.xfill$u, v = svd.xfill$v, cvg = cvg, iter = iter))
}

#' bls.cov
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param x (np)xN matrix of covariates
#' @param mu real number, direction of update for offset
#' @param alpha  direction of update for vector of regression parameters of length N
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param mu.tmp real number, current offset
#' @param alpha.tmp length N vector, current regression parameters
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param lambda2 positive number, regularization parameter for l1 norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @export
#' @import stats
#' @return A list with the following elements
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta a}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
#' @examples
#' y <- matrix(rnorm(6 * 10), nrow = 6)
#' x <- matrix(rnorm(60*2), nrow = 60)
#' mu <- 0
#' alpha <- rnorm(2)
#' theta <- matrix(rnorm(6 * 10), nrow = 6)
#' mut <- 0.1
#' alphat <- rnorm(2)
#' thetat <- matrix(rnorm(6 * 10), nrow = 6)
#' v <- rep("gaussian", 10)
#' t <- bls.cov(y, x, mu, alpha, theta, mut, alphat, thetat, lambda1 = 1, lambda2 = 1, var.type = v)
bls.cov <- function(y0, x, mu, alpha, theta, mu.tmp, alpha.tmp, theta.tmp,
                    b = 0.5, lambda1, lambda2, var.type, thresh = 1e-5, sc){
  #thresh = 0
  d <- dim(y0)
  n <- d[1]
  p <- d[2]
  q <- ncol(x)
  omega <- !is.na(y0)
  alpha.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha, nrow = n)
  alpha.tmp.mat <- matrix(matrix(as.numeric(x), nrow = n*p)%*%alpha.tmp, nrow = n)
  param <- mu + alpha.mat + theta
  param.tmp <- mu.tmp + alpha.tmp.mat + theta.tmp
  gaus.tmp <- (1 / 2) * sum(sc[, var.type == "gaussian"]*(y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"])^2,
                            na.rm = T)
  pois.tmp <- sum((- (y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
                                                 exp(param.tmp[, var.type == "poisson"])), na.rm = T)
  binom.tmp <- sum((- (y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
                                                 log(1 + exp(param.tmp[, var.type == "binary"]))), na.rm = T)
  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat.tmp <- sum((- (y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                                                    log(truc))), na.rm = T)
  } else cat.tmp <- 0
  d.tmp <- svd(theta.tmp)$d
  gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                        na.rm = T)
  pois <- sum((- (y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                                             exp(param[, var.type == "poisson"])), na.rm = T)
  binom <- sum((- (y0[, var.type == "binary"] * param[, var.type == "binary"]) +
                                             log(1 + exp(param[, var.type == "binary"]))), na.rm = T)
  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat <- sum((- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                log(truc))), na.rm = T)
  } else cat <- 0
  d <- svd(theta)$d
  t <- 1
  mu2 <- (1-t)*mu.tmp + t*mu
  alpha2 <- (1-t)*alpha.tmp + t*alpha
  theta2 <- (1-t)*theta.tmp + t*theta
  param2 <- (1-t)*param.tmp + t*param
  diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha))))
  number <- gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp) + sum(lambda2 * (t(abs(alpha.tmp))))
  while(diff < -abs(number)*thresh){
    t <- b*t
    mu2 <- (1-t)*mu.tmp + t*mu
    alpha2 <- (1-t)*alpha.tmp + t*alpha
    theta2 <- (1-t)*theta.tmp + t*theta
    param2 <- (1-t)*param.tmp + t*param
    gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param2[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum((- (y0[, var.type == "poisson"] * param2[, var.type == "poisson"]) +
                                               exp(param2[, var.type == "poisson"])), na.rm = T)
    binom <- sum((- (y0[, var.type == "binary"] * param2[, var.type == "binary"]) +
                                               log(1 + exp(param2[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if(sum(var.type=="categorical")>0){
      for(j in 1:sum(var.type=="categorical")){
        tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
        truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum((- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                  log(truc))), na.rm = T)
    } else cat <- 0
    d <- svd(theta2)$d
    diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha2))))
  }
  obj <- pois + gaus + binom + lambda1*d + sum(lambda2 * t(abs(alpha2)))
  return(list(mu = mu2, alpha = alpha2, theta = theta2, objective = obj, t=t))

}

#' bls.multi
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp matrix of observations
#' @param groups length n indicator factor of group memberships
#' @param mu real number new offset
#' @param alpha vector of length p*N(nb of groups) new regression parameter
#' @param theta nxp matrix new interaction matrix
#' @param mu.tmp real number current offset
#' @param alpha.tmp vector of length p*N(nb of groups) current regression parameter
#' @param theta.tmp nxp matrix current interaction matrix
#' @param b number in (0,1) factor by which step size is reduced
#' @param lambda1 positive number regularization parameter for nuclear norm penalty
#' @param lambda2 positive number regularization parameter for l1 norm penalty
#' @param var.type vector of length p indicating variable types (gaussian, binary, poisson)
#' @param thresh positive number congervence criterion
#'
#' @import stats
#' @export
#' @return A list with the following elements
#' \item{mu}{the offset}
#' \item{alpha}{a (nb groups) x (nb variables) matrix containing the group effects}
#' \item{theta a}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
#' @examples
#' y0 <- matrix(rnorm(6 * 10), nrow = 6)
#' m <- 0
#' a <- matrix(rnorm(30), nrow = 3)
#' t <- matrix(rnorm(6 * 10), nrow = 6)
#' mt <- 0.1
#' at <- matrix(rnorm(30), nrow = 3)
#' tt <- matrix(rnorm(6 * 10), nrow = 6)
#' v <- rep("gaussian", 10)
#' gps <- as.factor(c(1,1,2,2,3,3))
#' bls <- bls.multi(y0, gps, m, a, t, mt, at, tt, 1, 1, v)
bls.multi <- function(y0, groups, mu, alpha, theta, mu.tmp, alpha.tmp, theta.tmp,
                      lambda1, lambda2, var.type, thresh = 1e-5, b = 0.5){
  d <- dim(y0)
  n <- d[1]
  p <- d[2]
  groups <- as.factor(groups)
  N <- nlevels(groups)
  omega <- !is.na(y0)
  ncenters <- aggregate(rep(1, n), by = list(groups), sum)[,2]
  alpha.rep <- matrix(rep(as.matrix(alpha), rep(ncenters, p)), nrow = n)
  alpha.tmp.rep <- matrix(rep(as.matrix(alpha.tmp), rep(ncenters, p)), nrow = n)
  param <- mu + alpha.rep + theta
  param.tmp <- mu.tmp + alpha.tmp.rep + theta.tmp
  gaus.tmp <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"])^2,
                            na.rm = T)
  pois.tmp <- sum(- y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"] +
                                                 exp(param.tmp[, var.type == "poisson"]), na.rm = T)
  binom.tmp <- sum(- y0[, var.type == "binary"] * param.tmp[, var.type == "binary"] +
                                                 log(1 + exp(param.tmp[, var.type == "binary"])), na.rm = T)
  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat.tmp <- sum(- (y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                                                    log(truc)), na.rm = T)
  } else cat.tmp <- 0
  d.tmp <- svd(theta.tmp)$d
  gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                        na.rm = T)
  pois <- sum(-y0[, var.type == "poisson"]*param[, var.type == "poisson"] +
                                             exp(param[, var.type == "poisson"]), na.rm = T)
  binom <- sum(-y0[, var.type == "binary"]*param[, var.type == "binary"] +
                                             log(1 + exp(param[, var.type == "binary"])), na.rm = T)

  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat <- sum((- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                log(truc))), na.rm = T)
  } else cat <- 0
  d <- svd(theta)$d
  t <- 1
  mu2 <- (1-t)*mu.tmp + t*mu
  alpha2 <- (1-t)*alpha.tmp + t*alpha
  theta2 <- (1-t)*theta.tmp + t*theta
  param2 <- (1-t)*param.tmp + t*param
  diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha))))
  number <- gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp) + sum(lambda2 * (t(abs(alpha.tmp))))
  while(diff < -abs(number)*thresh/2){
    t <- b*t
    mu2 <- (1-t)*mu.tmp + t*mu
    alpha2 <- (1-t)*alpha.tmp + t*alpha
    theta2 <- (1-t)*theta.tmp + t*theta
    param2 <- (1-t)*param.tmp + t*param
    gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param2[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum((- (y0[, var.type == "poisson"] * param2[, var.type == "poisson"]) +
                                               exp(param2[, var.type == "poisson"])), na.rm = T)
    binom <- sum((- (y0[, var.type == "binary"] * param2[, var.type == "binary"]) +
                                               log(1 + exp(param2[, var.type == "binary"]))), na.rm = T)
    truc <- rep(0, n)
    if(sum(var.type=="categorical")>0){
      for(j in 1:sum(var.type=="categorical")){
        tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
        truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum((- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                  log(truc))), na.rm = T)
    } else cat <- 0
    d <- svd(theta2)$d

    diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d)) + sum(lambda2 * (t(abs(alpha.tmp)) - t(abs(alpha2))))
  }
  obj <- pois + gaus + binom + lambda1*d + sum(lambda2 * t(abs(alpha2)))
  return(list(mu = mu2, alpha = alpha2, theta = theta2, objective = obj, t=t))
}


#' bls
#' Performs backtracking line search along a pre-specified search direction
#' @param y0 nxp observations matrix
#' @param theta nxp matrix direction of update for matrix of interactions
#' @param theta.tmp nxp matrix, current matrix of interactions
#' @param b positive number in (0,1) factor by which the step size is reduced
#' @param lambda1 positive number, regularization parameter for nuclear norm penalty
#' @param var.type vector of length p indicating column types for y (gaussian, binary, poisson)
#' @param thresh positive number, convergence criterion
#' @export
#' @import stats
#' @return A list with the following elements
#' \item{theta}{(nb individuals) x (nb variables) matrix containing the individual effects}
#' \item{objective}{a vector containing the value of the objective function at every iteration}
#' \item{t}{the step size}
#' @examples
#' y <- matrix(rnorm(6 * 10), nrow = 6)
#' theta <- matrix(rnorm(6 * 10), nrow = 6)
#' thetat <- matrix(rnorm(6 * 10), nrow = 6)
#' v <- rep("gaussian", 10)
#' t <- bls.lr(y, theta, thetat, lambda1 = 1, var.type = v)
bls.lr <- function(y0, theta, theta.tmp, b = 0.5, lambda1, var.type, thresh = 1e-5){
  d <- dim(y0)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y0)
  param <- theta
  param.tmp <- theta.tmp
  gaus.tmp <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param.tmp[, var.type == "gaussian"])^2,
                            na.rm = T)
  pois.tmp <- sum(- (y0[, var.type == "poisson"] * param.tmp[, var.type == "poisson"]) +
                    exp(param.tmp[, var.type == "poisson"]), na.rm = T)
  binom.tmp <- sum(- (y0[, var.type == "binary"] * param.tmp[, var.type == "binary"]) +
                     log(1 + exp(param.tmp[, var.type == "binary"])), na.rm = T)
  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat.tmp <- sum(- (y0[, vt2 == "categorical"] * (param.tmp[, vt2 == "categorical"]) -
                                                log(truc)), na.rm = T)
  } else cat.tmp <- 0

  d.tmp <- svd(theta.tmp)$d
  gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param[, var.type == "gaussian"])^2,
                        na.rm = T)
  pois <- sum(- (y0[, var.type == "poisson"] * param[, var.type == "poisson"]) +
                exp(param[, var.type == "poisson"]), na.rm = T)
  binom <- sum(- (y0[, var.type == "binary"] * param[, var.type == "binary"]) +
                 log(1 + exp(param[, var.type == "binary"])), na.rm = T)

  truc <- rep(0, n)
  if(sum(var.type=="categorical")>0){
    for(j in 1:sum(var.type=="categorical")){
      tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat <- sum(- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                log(truc)), na.rm = T)
  } else cat <- 0
  d <- svd(theta)$d
  t <- 1
  theta2 <- (1-t)*theta.tmp + t*theta
  param2 <- (1-t)*param.tmp + t*param
  diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d))
  number <- gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp)
  while(diff < -abs(number)*thresh/2){
    t <- b*t
    theta2 <- (1-t)*theta.tmp + t*theta
    param2 <- (1-t)*param.tmp + t*param
    gaus <- (1 / 2) * sum((y0[, var.type == "gaussian"] - param2[, var.type == "gaussian"])^2,
                          na.rm = T)
    pois <- sum(- (y0[, var.type == "poisson"] * param2[, var.type == "poisson"]) +
                  exp(param2[, var.type == "poisson"]), na.rm = T)
    binom <- sum(- (y0[, var.type == "binary"] * param2[, var.type == "binary"]) +
                   log(1 + exp(param2[, var.type == "binary"])), na.rm = T)
    truc <- rep(0, n)
    if(sum(var.type=="categorical")>0){
      for(j in 1:sum(var.type=="categorical")){
        tt <- rowSums(exp(param[, which(var.type=="categorical")[j]:(which(var.type=="categorical")[j]+nlevel[which(var.type=="categorical")[j]]-1)]))
        truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type=="categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum(- (y0[, vt2 == "categorical"] * (param[, vt2 == "categorical"]) -
                                                  log(truc)), na.rm = T)
    } else cat <- 0
    d <- svd(theta2)$d
    diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + lambda1 * (sum(d.tmp) - sum(d))
  }
  obj <- pois + gaus + binom + lambda1*d
  return(list(theta = theta2, objective = obj, t=t))

}

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
covmat <- function(R,C,n,p){
  if(is.null(R)) {
    covs <- covmatC(C,n)
  } else if(is.null(C)) {
    covs <- covmatR(R,p)
  } else {
    R <- as.matrix(R)
    C <- as.matrix(C)
    dR <- dim(R)
    dC <- dim(C)
    K1 <- dR[2]
    K2 <- dC[2]
    covs <- cbind(do.call(rbind, replicate(nrow(C), R, simplify=FALSE)),
                  C[rep(seq_len(nrow(C)), each=nrow(R)),])
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
covmatR <- function(R, p){
  R <- as.matrix(R)
  dR <- dim(R)
  n <- dR[1]
  K1 <- dR[2]
  covs <- do.call(rbind, replicate(p, R, simplify=FALSE))
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
covmatC <- function(C, n){
  C <- as.matrix(C)
  dC <- dim(C)
  p <- dC[1]
  K2 <- dC[2]
  covs <- C[rep(seq_len(p), each=n),]
  return(covs)
}
