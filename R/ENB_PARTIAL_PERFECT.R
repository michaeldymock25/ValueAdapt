
#' @title enb_partial_perfect
#' @description Computes the expected net benefit of collecting partial perfect information using either a Monte-Carlo or non-parametric approximation method
#' @import mgcv
#' @import RhpcBLASctl
#' @import stats
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option, parameters and decision times t_1 and t_2
#' @param Theta_int Named matrix of parameter draws from prior/posterior distribution for parameters of interest
#' @param Theta_rem Named matrix of parameter draws from prior/posterior distribution for the remaining parameters (not of interest)
#' @param t Named vector of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H")
#' @param prop Vector containing current proportions of intervention use. Must sum to one.
#' @param cost Cost of sampling
#' @param method Approximation method. Either MC for Monte-Carlo, NP for non-parametric (default) or MM for moment matching. The moment matching method requires the evppi function to be run in advance using the non-parametric method to generate INB_partial.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param cond_args List of arguments to be passed to cond_fun(). Defaults to an empty list.
#' @param cond_fun A function that generates sets of Theta_rem conditional on the current values of Theta_int. The only argument should be a list containing the parameter draws and additional parameters. Only required for the Monte Carlo approximation method.
#' @param model Generalised additive regression model specification (formula). Only required for the non-parametric approximation method. Should be a function of the parameters of interest.
#' @return Expected net benefit of collecting partial perfect information
#' @examples
#' # two parameters (one parameter of interest), two decision options
#' D <- 2
#' U <- function(d, Theta_int, Theta_rem, t_1, t_2)
#'   sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int + (d == 2)*Theta_rem)
#' N <- 10000
#' Theta_int <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
#' Theta_rem <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' cond_fun <- function(args) rbeta(args$N, 2, 3)
#' enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "MC",
#'                     cond_fun = cond_fun)
#' enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "NP",
#'                     model = "s(theta_A)")
#'
#' # three parameters (two parameters of interest), three decision options
#' D <- 3
#' U <- function(d, Theta_int, Theta_rem, t_1, t_2)
#'   sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int[,"theta_A"] +
#'                            (d == 2)*Theta_int[,"theta_B"] +
#'                            (d == 3)*Theta_rem)
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' Theta_int <- Theta[,c("theta_A", "theta_B")]
#' Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' cond_fun <- function(args) rbeta(args$N, 2, 3)
#' enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "MC",
#'                     cond_fun = cond_fun)
#' enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "NP",
#'                     model = "te(theta_A, theta_B)")
#' @rdname enb_partial_perfect
#' @export
enb_partial_perfect <- function(D, U, Theta_int, Theta_rem, t, prop, cost, method = "NP", K = 10000,
                                cond_args = list(), cond_fun = NULL, model = NULL){

  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP")
  if(nrow(Theta_int) != nrow(Theta_rem)) stop("Theta_int and Theta_rem must contain the same number of draws")

  ## first compute the expected value of choosing now

  NB_now <- sapply(1:D, function(d) U(d, Theta_int, Theta_rem, t["C"] + 1, t["H"]))
  INB_now <- NB_now - NB_now[,1]
  value_now <- max(colMeans(INB_now))

  ## second compute the expected value during the trial

  if(sum(prop) != 1) stop("prop must sum to one")
  NB_during <- sapply(1:D, function(d) U(d, Theta_int, Theta_rem, t["C"] + 1, t["A"]))
  value_during <- mean(NB_during%*%prop - NB_during[,1])

  ## third compute the expected value of choosing after the trial

  N <- nrow(Theta_int)
  if(method == "MC"){
    if(N < K) stop("The number of parameter draws must be greater than or equal to K")
    Theta_int_redraw <- as.matrix(Theta_int[sample.int(N, size = K),])
    colnames(Theta_int_redraw) <- colnames(Theta_int)
    PPI <- sapply(1:K, function(k){
      cond_args_tmp <- c(cond_args, D = D, N = N, Theta_int_redraw[k,])
      Theta_rem_tmp <- cond_fun(cond_args_tmp)
      Theta_int_tmp <- do.call(rbind, lapply(1:N, function(n) Theta_int_redraw[k,]))
      NB_tmp <- sapply(1:D, function(d) U(d, Theta_int_tmp, Theta_rem_tmp, t["A"] + 1, t["H"]))
      INB_tmp <- NB_tmp - NB_tmp[,1]
      max(colMeans(INB_tmp))
    })
    value_after <- mean(PPI)

    ## finally compute the expected net benefit of collecting partial perfect information

    ENB_PARTIAL_PERFECT <- (value_during + value_after) - (value_now + cost)

    return(ENB_PARTIAL_PERFECT)

  } else if(method == "NP"){
    NB <- sapply(1:D, function(d) U(d, Theta_int, Theta_rem, t["A"] + 1, t["H"]))
    INB <- NB - NB[,1]
    INB_partial <- matrix(data = NA, nrow = N, ncol = D)
    INB_partial[,1] <- 0
    RhpcBLASctl::blas_set_num_threads(1)
    for(d in 2:D) INB_partial[,d] <- gam(update(formula(INB[, d] ~ .), formula(paste(".~", model))), data = as.data.frame(Theta_int))$fitted
    value_after <- mean(apply(INB_partial, 1, max))

    ## finally compute the expected net benefit of collecting partial perfect information

    ENB_PARTIAL_PERFECT <- (value_during + value_after) - (value_now + cost)

    return(list(ENB_PARTIAL_PERFECT = ENB_PARTIAL_PERFECT, INB_partial = INB_partial))
  }
}