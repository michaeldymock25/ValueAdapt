
#' @title enb_perfect
#' @description Computes the expected net benefit of collecting perfect information using a Monte-Carlo approximation method
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option, parameters and decision times t_1 and t_2
#' @param Theta Named matrix of parameter draws from prior/posterior distribution
#' @param t Named vector of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H")
#' @param prop Vector containing current proportions of intervention use. Must sum to one.
#' @param cost Cost of sampling
#' @return Expected net benefit of collecting perfect information
#' @examples
#' # one parameter, two decision options
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*(-1)^(d-1)*(Theta - 0.4)
#' N <- 10000
#' Theta <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta"))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' enb_perfect(D, U, Theta, t, prop, cost)
#'
#' # two parameters, two decision options
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 2, dimnames = list(NULL, c("theta_A", "theta_B")))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' enb_perfect(D, U, Theta, t, prop, cost)
#'
#' # three parameters, three decision options
#' D <- 3
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' enb_perfect(D, U, Theta, t, prop, cost)
#' @rdname enb_perfect
#' @export
enb_perfect <- function(D, U, Theta, t, prop, cost){

  ## first compute the expected value of choosing now

  NB_now <- sapply(1:D, function(d) U(d, Theta, t["C"] + 1, t["H"]))
  INB_now <- NB_now - NB_now[,1]
  value_now <- max(colMeans(INB_now))

  ## second compute the expected value during the trial

  if(sum(prop) != 1) stop("prop must sum to one")
  NB_during <- sapply(1:D, function(d) U(d, Theta, t["C"] + 1, t["A"]))
  value_during <- mean(NB_during%*%prop - NB_during[,1])

  ## third compute the expected value of choosing after the trial

  NB_after <- sapply(1:D, function(d) U(d, Theta, t["A"] + 1, t["H"]))
  INB_after <- NB_after - NB_after[,1]
  value_after <- mean(apply(INB_after, 1, max))

  ## finally compute the expected net benefit of collecting perfect information

  ENB_PERFECT <- (value_during + value_after - cost) - value_now

  return(ENB_PERFECT)
}