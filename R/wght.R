wght <- function(y, param, var.type) {
  if (var.type == "gaussian") {
    ytilde <- y
    vtilde2 <- 1
  } else if (var.type == "binomial") {
    vtilde2 <- 1 / 4
    ytilde <-
      y / vtilde2 - (exp(param) / (1 + exp(param))) / vtilde2 + param
  } else if (var.type == "poisson") {
    vtilde2 <- exp(param)
    ytilde <- (y - exp(param)) / vtilde2  + param
  }
  else {
    print(var.type)
    stop(
      "Incorrect type of variable. Should be 'gaussian', 'binomial', 'poisson' or 'categorical'."
    )
  }
  return(list(ytilde = ytilde, vtilde2 = vtilde2))
}
