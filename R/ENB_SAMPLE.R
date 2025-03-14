
#' @title enb_sample
#' @description Computes the expected net benefit of collecting sample information using either a Monte-Carlo, non-parametric or moment matching approximation method
#' @import EVSI
#' @import mgcv
#' @import parallel
#' @import stats
#' @import RhpcBLASctl
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option, parameters and decision times t_1 and t_2
#' @param Theta Named matrix of parameter draws from prior/posterior distribution
#' @param t Named vector of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H")
#' @param prop Vector containing current proportions of intervention use. Must sum to one.
#' @param cost Cost of sampling
#' @param method Approximation method. Either MC for Monte-Carlo, NP for non-parametric (default) or MM for moment matching. The moment matching method requires the evppi function to be run in advance using the non-parametric method to generate INB_partial. Must be MC if an underestimation correction is to be computed.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param samp_args List of arguments to be passed to samp_fun(). Defaults to an empty list.
#' @param samp_fun A function that generates a sample of data based on the parameter draws (posterior predictive draw). The only argument should be a list including parameter draws from Theta and any additional parameters. Returns a vector of data to be passed to post_fun().
#' @param post_args List of arguments to be passed to post_fun(). Defaults to an empty list.
#' @param post_fun A function that generates draws from the posterior distribution given the previous distribution and a sample of data generated by the samp_fun function. The only argument should be a list including a sample of data from samp_fun(), number of draws required from the posterior distribution (N) and any additional parameters. Only required for the Monte Carlo and moment matching approximation methods.
#' @param stat_fun A function to generate a low-dimensional summary statistic based on the sample of data generated by the samp_fun function. The only argument is a vector of sampled data. Only required for the non-parametric approximation method.
#' @param model Generalised additive regression model specification (formula). Only required for the non-parametric approximation method. Should be a function of the data exported by samp_fun().
#' @param INB_partial Samples of INB for the parameters of interest generated from the evppi function using the non-parametric approximation method. Only required for the moment matching approximation method.
#' @param Q Number of model reruns to estimate the expected variance of the posterior net benefit. Only required for the moment matching approximation method.
#' @param correct A named list containing information to compute a correction for the underestimation of the expected net benefit of sampling. Defaults to NULL to apply no correction.
#' @param num_cores Number of cores to utilise when generating posterior distributions. Only required for the Monte Carlo and moment matching approximation methods. Defaults to one.
#' @param num_threads Number of BLAS threads. Defaults to one.
#' @return Expected net benefit of collecting sample information
#' @examples
#' # one parameter, two decision options
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*(-1)^(d-1)*(Theta - 0.4)
#' N <- 10000
#' Theta <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta"))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' samp_fun <- function(args) c(x = rbinom(1, size = 10, prob = args$theta))
#' post_fun <- function(args) rbeta(args$N, 2 + args$x, 3 + 10 - args$x)
#' stat_fun <- function(x) x/10
#' enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
#' enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
#'            model = "s(x)")
#'
#' # two parameters (one parameter of interest), two decision options
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 2, dimnames = list(NULL, c("theta_A", "theta_B")))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' samp_fun <- function(args) c(x = rbinom(1, size = 10, prob = args$theta_A))
#' post_fun <- function(args) matrix(c(rbeta(args$N, 2 + args$x, 3 + 10 - args$x),
#'                                     rbeta(args$N, 2, 3)),
#'                                   nrow = args$N, ncol = 2,
#'                                   dimnames = list(NULL, c("theta_A", "theta_B")))
#' stat_fun <- function(x) x/10
#' enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
#' enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
#'            model = "s(x)")
#' U_enbppi <- function(d, Theta_int, Theta_rem, t_1, t_2)
#'   sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int + (d == 2)*Theta_rem)
#' Theta_int <- matrix(Theta[,"theta_A"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
#' Theta_rem <- matrix(Theta[,"theta_B"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
#' out <- enb_partial_perfect(D, U_enbppi, Theta_int, Theta_rem, t, prop, cost, method = "NP",
#'                            model = "s(theta_A)")
#' enb_sample(D, U, Theta, t, prop, cost, method = "MM", samp_fun = samp_fun, post_fun = post_fun,
#'            INB_partial = out$INB_partial)
#'
#' # three parameters (two parameters of interest), three decision options
#' D <- 3
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' N <- 10000
#' Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
#'                 nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' t <- c(C = 0, A = 1, H = 15)
#' prop <- rep(1/D, D)
#' cost <- 0
#' samp_fun <- function(args) c(x_A = rbinom(1, size = 10, prob = args$theta_A),
#'                              x_B = rbinom(1, size = 10, prob = args$theta_B))
#' post_fun <- function(args) matrix(c(rbeta(args$N, 2 + args$x_A, 3 + 10 - args$x_A),
#'                                     rbeta(args$N, 2 + args$x_B, 3 + 10 - args$x_B),
#'                                     rbeta(args$N, 2, 3)),
#'                                   nrow = args$N, ncol = 3,
#'                                   dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
#' stat_fun <- function(x) x/10
#' enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
#' enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
#'            model = "te(x_A, x_B)")
#' U_enbppi <- function(d, Theta_int, Theta_rem, t_1, t_2)
#'   sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int[,"theta_A"] +
#'                            (d == 2)*Theta_int[,"theta_B"] +
#'                            (d == 3)*Theta_rem)
#' Theta_int <- Theta[,c("theta_A", "theta_B")]
#' Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
#' out <- enb_partial_perfect(D, U_enbppi, Theta_int, Theta_rem, t, prop, cost, method = "NP",
#'                            cond_fun = cond_fun, model = "te(theta_A, theta_B)")
#' enb_sample(D, U, Theta, t, prop, cost, method = "MM", samp_fun = samp_fun,
#'            post_fun = post_fun, INB_partial = out$INB_partial)
#' @rdname enb_sample
#' @export
enb_sample <- function(D, U, Theta, t, prop, cost, method = "NP", K = 10000,
                       samp_args = list(), samp_fun = NULL, post_args = list(), post_fun = NULL,
                       stat_fun = NULL, model = NULL, INB_partial = NULL, Q = 50,
                       correct = NULL, num_cores = 1, num_threads = 1){

  RhpcBLASctl::blas_set_num_threads(num_threads)
  if(!(method %in% c("MC", "NP", "MM"))) stop("Method must be specified as MC, NP or MM")
  if(!is.null(correct) & method != "MC") stop("Method must be MC to calculate the correction term")

  ## first compute the expected value of choosing now

  NB_now <- sapply(1:D, function(d) U(d, Theta, t["C"] + 1, t["H"]))
  INB_now <- NB_now - NB_now[,1]
  value_now <- max(colMeans(INB_now))

  ## second compute the expected value during the trial

  if(sum(prop) != 1) stop("prop must sum to one")
  NB_during <- sapply(1:D, function(d) U(d, Theta, t["C"] + 1, t["A"]))
  value_during <- mean(NB_during%*%prop - NB_during[,1])

  ## third compute the expected value of choosing after the trial

  N <- nrow(Theta)
  if(method == "MC"){
    if(N < K) stop("The number of parameter draws must be greater than or equal to K")
    Theta_redraw <- as.matrix(Theta[sample.int(N, size = K),])
    colnames(Theta_redraw) <- colnames(Theta)
    samp_out <- lapply(1:K, function(k){
                   samp_args_tmp <- c(samp_args, Theta_redraw[k,])
                   samp_fun(samp_args_tmp)
    })
    Theta_tmp <- mclapply(1:K,
                          function(k){
                            post_args_tmp <- c(post_args, D = D, N = N, samp_out[[k]])
                            post_fun(post_args_tmp)
                          },
                          mc.cores = num_cores)
    SI <- sapply(1:K, function(k){
      NB_tmp <- sapply(1:D, function(d) U(d, Theta_tmp[[k]], t["A"] + 1, t["H"]))
      INB_tmp <- NB_tmp - NB_tmp[,1]
      max(colMeans(INB_tmp))
    })
    value_after <- mean(SI)
  } else if(method == "NP"){
    NB <- sapply(1:D, function(d) U(d, Theta, t["A"] + 1, t["H"]))
    INB <- NB - NB[,1]
    samp_out <- lapply(1:N, function(n){
                   Theta_tmp <- Theta[n,]
                   names(Theta_tmp) <- colnames(Theta)
                   samp_args_tmp <- c(samp_args, Theta_tmp)
                   samp_fun(samp_args_tmp)
    })
    summ_stats <- t(matrix(sapply(samp_out, stat_fun), ncol = N))
    colnames(summ_stats) <- names(samp_out[[1]])
    g_hat <- matrix(data = NA, nrow = N, ncol = D)
    g_hat[,1] <- 0
    for(d in 2:D) g_hat[,d] <- gam(update(formula(INB[, d] ~ .), formula(paste(".~", model))), data = as.data.frame(summ_stats))$fitted
    value_after <- mean(apply(g_hat, 1, max))
  } else if(method == "MM"){
    NB <- sapply(1:D, function(d) U(d, Theta, t["A"] + 1, t["H"]))
    INB <- NB - NB[,1]
    Theta_redraw <- as.matrix(gen.quantiles(parameter = colnames(Theta), param.mat = Theta, Q = Q))
    samp_out <- lapply(1:Q, function(q){
                   samp_args_tmp <- c(samp_args, Theta_redraw[q,])
                   samp_fun(samp_args_tmp)
    })
    Theta_tmp <- mclapply(1:Q,
                          function(q){
                            post_args_tmp <- c(post_args, D = D, N = N, samp_out[[q]])
                            post_fun(post_args_tmp)
                          },
                          mc.cores = num_cores)
    var_est <- lapply(1:Q, function(q){
      NB_tmp <- sapply(1:D, function(d) U(d, Theta_tmp[[q]], t["A"] + 1, t["H"]))
      INB_tmp <- NB_tmp - NB_tmp[,1]
      var(INB_tmp[,-1])
    })
    prior_var <- var(as.matrix(INB[,-1]))
    mu_mn <- colMeans(as.matrix(INB[,-1]))
    mu_var <- prior_var - Reduce("+", var_est)/Q
    if(D <= 2){
      INB_rescaled <- as.matrix((INB_partial[,-1] - mu_mn) * 1/sd(INB_partial[,-1]) * as.numeric(sqrt(mu_var)) + mu_mn)
    } else {
      mu_var_eigen <- eigen(mu_var)
      mu_var_sqrt <- mu_var_eigen$vectors %*% diag(sqrt(mu_var_eigen$values)) %*% t(mu_var_eigen$vectors)
      prior_var_eigen <- eigen(prior_var)
      prior_var_sqrt_inv <- chol2inv(chol(prior_var_eigen$vectors %*% diag(sqrt(prior_var_eigen$values)) %*% t(prior_var_eigen$vectors)))

      INB_rescaled <- t(t(t(t(INB_partial[,-1]) - mu_mn) %*% prior_var_sqrt_inv %*% mu_var_sqrt) + mu_mn)
    }
    value_after <- mean(apply(INB_rescaled, 1, function(x) max(0, x)))
  }

  ## fourth compute the underestimation correction term if specified

  if(!is.null(correct)){
    if(is.null(correct[["n_analyses_rem"]])) correct[["n_analyses_rem"]] <- correct[["n_analyses"]] - 1
    if(correct[["n_analyses_rem"]] > 0){
      correct[["n_analyses_rem"]] <- correct[["n_analyses_rem"]] - 1
      analysis <- correct[["n_analyses"]] - correct[["n_analyses_rem"]]
      if("t_update" %in% names(correct)){
        if(!("t_update_args" %in% names(correct))) correct[["t_update_args"]] <- list()
        correct[["t_update_args"]][["analysis"]] <- analysis
        t <- correct[["t_update"]](correct[["t_update_args"]])
      }
      if("cost_update" %in% names(correct)){
        if(!("cost_update_args" %in% names(correct))) correct[["cost_update_args"]] <- list()
        correct[["cost_update_args"]][["analysis"]] <- analysis
        cost <- correct[["cost_update"]](correct[["cost_update_args"]])
      }
      if("samp_args_update" %in% names(correct)){
        if(!("samp_args_update_args" %in% names(correct))) correct[["samp_args_update_args"]] <- list()
        correct[["samp_args_update_args"]][["analysis"]] <- analysis
        samp_args <- correct[["samp_args_update"]](correct[["samp_args_update_args"]])
      }
      correction_term_vec <- sapply(1:K, function(k){
                                if(!("post_args_update_args" %in% names(correct))) correct[["post_args_update_args"]] <- list()
                                for(nm in names(samp_out[[k]])) correct[["post_args_update_args"]][[nm]] <- samp_out[[k]][[nm]]
                                for(nm in names(post_args)) correct[["post_args_update_args"]][[nm]] <- post_args[[nm]]
                                correct[["post_args_update_args"]][["analysis"]] <- analysis
                                post_args <- correct[["post_args_update"]](correct[["post_args_update_args"]])
                                enb_sample(D = D, U = U, Theta = Theta_tmp[[k]], t = t, prop = prop, cost = cost, method = method, K = K,
                                           samp_args = samp_args, samp_fun = samp_fun, post_args = post_args, post_fun = post_fun,
                                           correct = correct)})
      correction_term <- mean(pmax(correction_term_vec, 0))
    } else {
      correction_term <- 0
    }
  } else {
    correction_term <- 0
  }

  ## finally compute the expected net benefit of collecting sample information

  ENB_SAMPLE <- (value_during + value_after + correction_term) - (value_now + cost)

  return(ENB_SAMPLE)
}