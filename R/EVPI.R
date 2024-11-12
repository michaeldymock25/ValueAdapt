
#' @title evpi
#' @description Computes the expected value of perfect information using a Monte-Carlo approximation method
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option and parameters
#' @param Theta Named matrix of parameter draws from prior/posterior distribution
#' @return Expected value of perfect information
#' @examples
#' # one parameter, two decision options
#' D <- 2
#' U <- function(d, Theta) (-1)^(d-1)*(Theta - 0.4)
#' N <- 10000
#' Theta <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta"))
#' evpi(D, U, Theta)
#'
#' # two parameters, two decision options
#' D <- 2
#' U <- function(d, Theta) Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 2, dimnames = list(NULL, c("theta_A", "theta_B")))
#' evpi(D, U, Theta)
#'
#' # three parameters, three decision options
#' D <- 3
#' U <- function(d, Theta) Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' evpi(D, U, Theta)
#' @rdname evpi
#' @export
evpi <- function(D, U, Theta){

  ## compute incremental net benefit using utility function U for each decision option d

  NB <- sapply(1:D, function(d) U(d, Theta))
  INB <- NB - NB[,1]

  ## estimate expected value of perfect information

  EVPI <- mean(apply(INB, 1, max)) - max(colMeans(INB))

  return(EVPI)
}