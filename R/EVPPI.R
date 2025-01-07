
#' @title evppi
#' @description Computes the expected value of partial perfect information using either a Monte-Carlo or non-parametric approximation method
#' @import mgcv
#' @import RhpcBLASctl
#' @import stats
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option and parameters
#' @param Theta_int Named matrix of parameter draws from prior/posterior distribution for parameters of interest
#' @param Theta_rem Named matrix of parameter draws from prior/posterior distribution for the remaining parameters (not of interest)
#' @param method Approximation method. Either MC for Monte-Carlo or NP for non-parametric (default)
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param cond_args List of arguments to be passed to cond_fun(). Defaults to an empty list.
#' @param cond_fun A function that generates sets of Theta_rem conditional on the current values of Theta_int. The only argument should be a list containing the parameter draws and additional parameters. Only required for the Monte Carlo approximation method.
#' @param model Generalised additive regression model specification (formula). Only required for the non-parametric approximation method. Should be a function of the parameters of interest.
#' @return Expected value of partial perfect information. If using the non-parametric method, a list will be returned containing the expected value of partial perfect information in addition to the partial incremental net benefit.
#' @examples
#' # two parameters (one parameter of interest), two decision options
#' D <- 2
#' U <- function(d, Theta_int, Theta_rem) (d == 1)*Theta_int + (d == 2)*Theta_rem
#' N <- 10000
#' Theta_int <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
#' Theta_rem <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
#' cond_fun <- function(args) rbeta(args$N, 2, 3)
#' evppi(D, U, Theta_int, Theta_rem, method = "MC", cond_fun = cond_fun)
#' evppi(D, U, Theta_int, Theta_rem, method = "NP", model = "s(theta_A)")
#'
#' # three parameters (two parameters of interest), three decision options
#' D <- 3
#' U <- function(d, Theta_int, Theta_rem) (d == 1)*Theta_int[,"theta_A"] +
#'                                        (d == 2)*Theta_int[,"theta_B"] +
#'                                        (d == 3)*Theta_rem
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' Theta_int <- Theta[,c("theta_A", "theta_B")]
#' Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
#' cond_fun <- function(args) rbeta(args$N, 2, 3)
#' evppi(D, U, Theta_int, Theta_rem, method = "MC", cond_fun = cond_fun)
#' evppi(D, U, Theta_int, Theta_rem, method = "NP", model = "te(theta_A, theta_B)")
#' @rdname evppi
#' @export
evppi <- function(D, U, Theta_int, Theta_rem, method = "NP", K = 10000, cond_args = list(), cond_fun = NULL, model = NULL){

  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP")
  if(nrow(Theta_int) != nrow(Theta_rem)) stop("Theta_int and Theta_rem must contain the same number of draws")

  ## compute incremental net benefit using utility function U for each decision option d

  NB <- sapply(1:D, function(d) U(d, Theta_int, Theta_rem))
  INB <- NB - NB[,1]

  ## estimate the expected value of partial perfect information

  N <- nrow(Theta_int)
  if(method == "MC"){
    if(N < K) stop("The number of parameter draws must be greater than or equal to K")
    Theta_int_redraw <- as.matrix(Theta_int[sample.int(N, size = K),])
    colnames(Theta_int_redraw) <- colnames(Theta_int)
    PPI <- sapply(1:K, function(k){
      cond_args_tmp <- c(cond_args, D = D, N = N, Theta_int_redraw[k,])
      Theta_rem_tmp <- cond_fun(cond_args_tmp)
      Theta_int_tmp <- do.call(rbind, lapply(1:N, function(n) Theta_int_redraw[k,]))
      NB_tmp <- sapply(1:D, function(d) U(d, Theta_int_tmp, Theta_rem_tmp))
      INB_tmp <- NB_tmp - NB_tmp[,1]
      max(colMeans(INB_tmp))
    })
    EVPPI <- mean(PPI) - max(colMeans(INB))

    return(EVPPI)

  } else if(method == "NP"){
    INB_partial <- matrix(data = NA, nrow = N, ncol = D)
    INB_partial[,1] <- 0
    RhpcBLASctl::blas_set_num_threads(1)
    for(d in 2:D) INB_partial[,d] <- gam(update(formula(INB[, d] ~ .), formula(paste(".~", model))), data = as.data.frame(Theta_int))$fitted
    EVPPI <- mean(apply(INB_partial, 1, max)) - max(colMeans(INB_partial))

    return(list(EVPPI = EVPPI, INB_partial = INB_partial))
  }
}