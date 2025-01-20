
#' @title sim_betabinomial_trial
#' @import parallel
#' @description Simulates a series of sequential trials using the value-driven adaptive design for a beta-binomial model
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option and parameters
#' @param Theta Vector of true parameter values to use for simulation. Must be of length D.
#' @param n_analyses Maximum number of trial analyses (must be at least two)
#' @param t Named matrix of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H") for each analysis. Columns are times and rows are analyses (with n_analyses rows).
#' @param prop Vector containing current proportions of intervention use. Must sum to one.
#' @param cost Vector of costs to proceed between analyses (not including trial start-up). Must be of length n_analyses.
#' @param prior_par Matrix of parameters to specify the beta prior distribution with each row representing a different parameter contained in Theta.
#' @param n_samp Number of observations to sample between analyses for each decision option
#' @param n_sims Number of simulated trials to generate. Defaults to one.
#' @param n_draws Number of parameter samples to draw from prior/posterior distributions. Defaults to 10,000.
#' @param method Approximation method. Either MC for Monte-Carlo or NP for non-parametric (default). The moment matching method is not currently supported. Must be MC if an underestimation correction is to be computed.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param apply_correct If TRUE, a correction for the underestimation of the expected net benefit of sampling is applied. Defaults to FALSE to apply no correction.
#' @param n_cores Number of cores to use if running simulations in parallel. Note that the num_cores and threads arguments of enb_sample() is set to one. Deafults to one.
#' @return Matrix containing the expected value of sample information computed at each analysis for each simulated trial
#' @examples
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' Theta <- c(0.5, 0.6)
#' n_analyses <- 3
#' t <- matrix(c(0, 1, 15,
#'               1, 2, 15,
#'               2, 3, 15),
#'             nrow = n_analyses, ncol = 3, byrow = TRUE,
#'             dimnames = list(NULL, c("C", "A", "H")))
#' prop <- rep(1/D, D)
#' cost <- rep(0, n_analyses)
#' prior_par <- matrix(rep(c(2, 3), D), nrow = D, byrow = TRUE)
#' sim_betabinomial_trial(D, U, Theta, n_analyses, t, prop, cost,
#'                        prior_par, n_samp = 50, n_sims = 4)
#' @rdname sim_betabinomial_trial
#' @export
sim_betabinomial_trial <- function(D, U, Theta, n_analyses, t, prop, cost, prior_par, n_samp, n_sims = 1,
                                   n_draws = 10000, method = "NP", K = 10000, apply_correct = FALSE, n_cores = 1){

  if(length(Theta) != D) stop("Theta must be a vector of length D")
  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP. Note that the moment matching method is not currently supported")
  if(nrow(t) != n_analyses) stop("Matrix t must have rows for the number of analyses")
  if(length(cost) != n_analyses) stop("Cost vector must be of length the number of analyses")
  if(apply_correct & method != "MC") stop("Method must be MC to calculate the correction term")

  dat_trial <- rbinom(n_sims*n_analyses*D, size = n_samp, Theta)
  dat_trial <- array(dat_trial,
                     dim = c(D, n_analyses, n_sims),
                     dimnames = list("Decision" = paste0("d", 1:D), "Analysis" = 1:n_analyses, "Simulation" = 1:n_sims))
  dat_trial <- apply(dat_trial, c("Decision", "Simulation"), cumsum)

  prior <- mclapply(1:n_sims,
                    function(sim) lapply(1:D, function(d) rbeta(n_draws, prior_par[d,1], prior_par[d,2])),
                    mc.cores = n_cores)
  prior <- array(unlist(prior),
                 dim = c(n_draws, D, n_sims),
                 dimnames = list("Draw" = 1:n_draws, "Decision" = paste0("p_", 1:D), "Simulation" = 1:n_sims))

  posterior <- mclapply(1:n_sims,
                        function(sim){
                          lapply(1:n_analyses, function(analysis) lapply(1:D, function(d)
                            rbeta(n_draws,
                                  prior_par[d,1] + dat_trial[analysis, d, sim],
                                  prior_par[d,2] + n_samp*analysis - dat_trial[analysis, d, sim])))
                        },
                        mc.cores = n_cores)
  posterior <- array(unlist(posterior),
                     dim = c(n_draws, D, n_analyses, n_sims),
                     dimnames = list("Draw" = 1:n_draws, "Decision" = paste0("p_", 1:D), "Analysis" = 1:n_analyses, "Simulation" = 1:n_sims))

  samp_args <- list(n = n_samp)
  samp_fun <- function(args){
    out <- c(args[["n"]], sapply(1:D, function(d) rbinom(1, size = args[["n"]], prob = args[[paste0("p_", d)]])))
    names(out) <- c("n", paste0("y_", 1:D))
    out
  }

  post_args <- c(list(prior_par), rep(list(0), D + 1))
  names(post_args) <- c("prior_par", "n_prev", paste0("y_", 1:D, "_prev"))
  post_fun <- function(args){
    matrix(unlist(lapply(1:D, function(d)
                     rbeta(args[["N"]],
                           args[["prior_par"]][d,1] + (args[[paste0("y_", d, "_prev")]] + args[[paste0("y_", d)]]),
                           args$prior_par[d,2] + (args[["n_prev"]] + args[["n"]]) -
                             (args[[paste0("y_", d, "_prev")]] + args[[paste0("y_", d)]])))),
           nrow = args$N, ncol = D, dimnames = list(NULL, paste0("p_", 1:D)))
  }

  if(apply_correct){
    correct <- list(n_analyses = n_analyses,
                    t_update_args = list(t = t),
                    t_update = function(t_update_args) t_update_args[["t"]][t_update_args[["analysis"]],],
                    cost_update_args = list(cost = cost),
                    cost_update = function(cost_update_args) cost_update_args[["cost"]][cost_update_args[["analysis"]]],
                    post_args_update = function(post_args_update_args){
                                         out <- c(list(post_args_update_args$prior_par,
                                                       post_args_update_args$n_prev + post_args_update_args$n),
                                                  lapply(1:D, function(d) post_args_update_args[[paste0("y_", d, "_prev")]] +
                                                                          post_args_update_args[[paste0("y_", d)]]))
                                         names(out) <- c("prior_par", "n_prev", paste0("y_", 1:D, "_prev"))
                                         out})
  } else {
    correct <- NULL
  }

  ENB_SAMPLE <- mclapply(1:n_sims,
                         function(sim){
                           lapply(1:n_analyses, function(analysis){
                              if(analysis == 1){
                                enb_sample(D = D, U = U, Theta = prior[,,sim], t = t[analysis,],
                                           prop = prop, cost = cost[analysis], method = method, K = K,
                                           samp_args = samp_args, samp_fun = samp_fun,
                                           post_args = post_args, post_fun = post_fun,
                                           stat_fun = function(x) x, model = paste0("s(", paste0("y_", 1:D), ")", collapse = " + "),
                                           correct = correct)
                              } else {
                                enb_sample(D = D, U = U, Theta = posterior[,,analysis-1,sim], t = t[analysis,],
                                           prop = prop, cost = cost[analysis], method = method, K = K,
                                           samp_args = samp_args, samp_fun = samp_fun,
                                           post_args = post_args, post_fun = post_fun,
                                           stat_fun = function(x) x, model = paste0("s(", paste0("y_", 1:D), ")", collapse = " + "),
                                           correct = correct)
                              }
                           })
                         },
                         mc.cores = n_cores)

  ENB_SAMPLE <- array(unlist(ENB_SAMPLE),
                      dim = c(n_analyses, n_sims),
                      dimnames = list("Analysis" = 1:n_analyses, "Simulation" = 1:n_sims))

  return(ENB_SAMPLE)
}