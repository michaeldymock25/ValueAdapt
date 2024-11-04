
#' @title evppi
#' @description Computes the expected value of partial perfect information using either a Monte-Carlo or non-parametric approximation method
#' @import mgcv
#' @import stats
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option and parameters
#' @param Theta_int Named matrix of parameter draws from prior/posterior distribution for parameters of interest
#' @param Theta_rem Named matrix of parameter draws from prior/posterior distribution for the remaining parameters (not of interest)
#' @param method Approximation method. Either MC for Monte-Carlo or NP for non-parametric (default)
#' @param J Number of inner Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param cond_fun A function that generates sets of Theta_rem conditional on the current values of Theta_int. The first argument must be J and the second must be a parameter vector. Only required for the Monte Carlo approximation method.
#' @param model Generalised additive regression model specification (formula). Only required for the non-parametric approximation method.
#' @return Expected value of partial perfect information. If using the non-parametric method, a list will be returned containing the expected value of partial perfect information in addition to the partial incremental net benefit.
#' @examples
#' # two parameters (one parameter of interest), two decision options
#' D <- 2
#' U <- function(d, Theta_int, Theta_rem) (d == 1)*Theta_int +
#'                                        (d == 2)*Theta_rem
#' N <- 10000
#' Theta_int <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
#' Theta_rem <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
#' cond_fun <- function(J, Theta_int) rbeta(J, 2, 3)
#' evppi(D, U, Theta_int, Theta_rem, method = "MC", cond_fun = cond_fun)
#' evppi(D, U, Theta_int, Theta_rem, method = "NP", model = "s(theta_A)")
#'
#' # three parameters (two parameters of interest), three decision options
#' D <- 3
#' U <- function(d, Theta_int, Theta_rem) (d == 1)*Theta_int[,"theta_A"] +
#'                                        (d == 2)*Theta_int[,"theta_B"] +
#'                                        (d == 3)*Theta_rem
#' N <- 10000
#' Theta_int <- Theta[,c("theta_A", "theta_B")]
#' Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
#' cond_fun <- function(J, Theta_int) rbeta(J, 2, 3)
#' evppi(D, U, Theta_int, Theta_rem, method = "MC", cond_fun = cond_fun)
#' evppi(D, U, Theta_int, Theta_rem, method = "NP", model = "te(theta_A, theta_B)")
#' @rdname evppi
#' @export
evppi <- function(D, U, Theta_int, Theta_rem, method = "NP", J = 10000, K = 10000, cond_fun = NULL, model = NULL){

  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP")
  if(nrow(Theta_int) != nrow(Theta_rem)) stop("Theta_int and Theta_rem must contain the same number of draws")
  N <- nrow(Theta_int)

  ## compute incremental net benefit using utility function U for each decision option d

  NB <- sapply(1:D, function(d) U(d, Theta_int, Theta_rem))
  INB <- NB - NB[,1]

  ## estimate the expected value of partial perfect information

  if(method == "MC"){
    if(N < K) stop("The number of parameter draws must be greater than or equal to K")
    Theta_int_redraw <- as.matrix(Theta_int[sample.int(N, size = K),])
    PPI <- sapply(1:K, function(k){
      Theta_rem_tmp <- cond_fun(J, Theta_int_redraw[k,])
      Theta_int_tmp <- do.call(rbind, lapply(1:J, function(j) Theta_int_redraw[k,]))
      NB_tmp <- sapply(1:D, function(d) U(d, Theta_int_tmp, Theta_rem_tmp))
      INB_tmp <- NB_tmp - NB_tmp[,1]
      max(colMeans(INB_tmp))
    })
    EVPPI <- mean(PPI) - max(colMeans(INB))

    return(EVPPI)

  } else if(method == "NP"){
    INB_partial <- matrix(data = NA, nrow = N, ncol = D)
    INB_partial[,1] <- 0
    for(d in 2:D) INB_partial[,d] <- gam(update(formula(INB[, d] ~ .), formula(paste(".~", model))), data = as.data.frame(Theta_int))$fitted
    EVPPI <- mean(apply(INB_partial, 1, max)) - max(colMeans(INB_partial))

    return(list(EVPPI = EVPPI, INB_partial = INB_partial))
  }
}