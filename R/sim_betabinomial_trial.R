
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
#' @param n_samp Number of observations to sample between analyses for each decision option between each analysis. Currently only supports a scalar input.
#' @param n_sims Number of simulated trials to generate. Defaults to one.
#' @param n_draws Number of parameter samples to draw from prior/posterior distributions. Defaults to 10,000.
#' @param method Approximation method. Either MC for Monte-Carlo or NP for non-parametric (default). The moment matching method is not currently supported. Must be MC if an underestimation correction is to be computed.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param apply_correct If TRUE, a correction for the underestimation of the expected net benefit of sampling is applied. Defaults to FALSE to apply no correction.
#' @param n_cores Number of cores to use if running simulations in parallel. Note that the num_cores and threads arguments of enb_sample() is set to one. Deafults to one.
#' @return List containing arrays of the chosen decision options (DEC), expected net benefit of perfect information (ENBP), expected net benefit of sample information (ENBS) and run time (RUN_TIME) computed at each analysis for each simulated trial
#' @examples
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(-(t_1:(t_2-1))))*Theta[,d]
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
#' n_samp <- 50
#' n_sims <- 4
#' sim_betabinomial_trial(D, U, Theta, n_analyses, t, prop, cost,
#'                        prior_par, n_samp, n_sims)
#' @rdname sim_betabinomial_trial
#' @export
sim_betabinomial_trial <- function(D, U, Theta, n_analyses, t, prop, cost, prior_par, n_samp, n_sims = 1,
                                   n_draws = 10000, method = "NP", K = 10000, apply_correct = FALSE, n_cores = 1){

  if(length(Theta) != D) stop("Theta must be a vector of length D")
  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP. Note that the moment matching method is not currently supported")
  if(nrow(t) != n_analyses) stop("Matrix t must have the same number of rows as the number of analyses")
  if(length(cost) != n_analyses) stop("Cost vector must be the same length as the number of analyses")
  if(apply_correct & method != "MC") stop("Method must be MC to calculate the correction term")

  dat_trial <- as.vector(sapply(1:n_sims, function(sim)
                            sapply(0:(n_analyses - 1), function(analysis){
                              if(analysis == 0){
                                rep(0, D)
                              } else {
                                rbinom(D, size = n_samp, Theta)
                              }
                            })))
  dat_trial <- array(dat_trial,
                     dim = c(D, n_analyses, n_sims),
                     dimnames = list("Decision" = paste0("d", 1:D), "Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))
  dat_trial <- apply(dat_trial, c("Decision", "Simulation"), cumsum)

  psa <- mclapply(1:n_sims,
                  function(sim){
                    lapply(1:n_analyses, function(analysis)
                       lapply(1:D, function(d)
                           rbeta(n_draws,
                                 prior_par[d,1] + dat_trial[analysis, d, sim],
                                 prior_par[d,2] + n_samp*(analysis - 1) - dat_trial[analysis, d, sim])))
                  },
                  mc.cores = n_cores)
  psa <- array(unlist(psa),
               dim = c(n_draws, D, n_analyses, n_sims),
               dimnames = list("Draw" = 1:n_draws, "Decision" = paste0("p", 1:D), "Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))

  samp_args <- list(n = n_samp)
  samp_fun <- function(args){
    out <- c(args[["n"]], sapply(1:D, function(d) rbinom(1, size = args[["n"]], prob = args[[paste0("p", d)]])))
    names(out) <- c("n", paste0("d", 1:D))
    out
  }

  post_args <- c(list(prior_par), rep(list(0), D + 1))
  names(post_args) <- c("prior_par", "n_prev", paste0("d", 1:D, "_prev"))
  post_fun <- function(args){
    matrix(unlist(lapply(1:args[["D"]], function(d)
                     rbeta(args[["N"]],
                           args[["prior_par"]][d,1] + (args[[paste0("d", d, "_prev")]] + args[[paste0("d", d)]]),
                           args$prior_par[d,2] + (args[["n_prev"]] + args[["n"]]) - (args[[paste0("d", d, "_prev")]] + args[[paste0("d", d)]])))),
           nrow = args[["N"]], ncol = args[["D"]], dimnames = list(NULL, paste0("p", 1:args[["D"]])))
  }

  ENB <- mclapply(1:n_sims,
                  function(sim){
                    lapply(1:n_analyses, function(analysis){
                       if(apply_correct){
                         correct <- list(n_analyses = n_analyses - (analysis - 1),
                                         t_update_args = list(t = t[analysis:n_analyses,]),
                                         t_update = function(t_update_args) t_update_args[["t"]][t_update_args[["analysis"]],],
                                         cost_update_args = list(cost = cost[analysis:n_analyses]),
                                         cost_update = function(cost_update_args) cost_update_args[["cost"]][cost_update_args[["analysis"]]],
                                         post_args_update = function(post_args_update_args){
                                           out <- c(list(post_args_update_args[["prior_par"]],
                                                         post_args_update_args[["n_prev"]] + post_args_update_args[["n"]]),
                                                    lapply(1:D, function(d)
                                                       post_args_update_args[[paste0("d", d, "_prev")]] +
                                                       post_args_update_args[[paste0("d", d)]]))
                                           names(out) <- c("prior_par", "n_prev", paste0("d", 1:D, "_prev"))
                                           out})
                       } else {
                         correct <- NULL
                       }
                       Theta_tmp <- psa[,,analysis,sim]
                       t_tmp <- t[analysis,]
                       dec <- which.max(colMeans(sapply(1:D, function(d) U(d, Theta_tmp, t_tmp["C"], t_tmp["H"]))))
                       enbp <- enb_perfect(D = D, U = U, Theta = Theta_tmp, t = t_tmp, prop = prop, cost = cost[analysis])
                       start_time <- Sys.time()
                       enbs <- enb_sample(D = D, U = U, Theta = Theta_tmp, t = t_tmp, prop = prop, cost = cost[analysis], method = method,
                                          K = K, samp_args = samp_args, samp_fun = samp_fun, post_args = post_args, post_fun = post_fun,
                                          stat_fun = function(x) x, model = paste0("s(", paste0("d", 1:D), ")", collapse = " + "),
                                          correct = correct)
                       end_time <- Sys.time()
                       run_time <- difftime(end_time, start_time, units = "mins")
                       return(list(dec = dec, enbp = enbp, enbs = enbs, run_time = run_time))
                    })
                  },
                  mc.cores = n_cores)

  DEC <- array(sapply(ENB, function(analysis) sapply(analysis, function(x) x$dec)),
                      dim = c(n_analyses, n_sims),
                      dimnames = list("Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))

  ENB_PERFECT <- array(sapply(ENB, function(analysis) sapply(analysis, function(x) x$enbp)),
                       dim = c(n_analyses, n_sims),
                       dimnames = list("Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))

  ENB_SAMPLE <- array(sapply(ENB, function(analysis) sapply(analysis, function(x) x$enbs)),
                      dim = c(n_analyses, n_sims),
                      dimnames = list("Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))

  RUN_TIME <- array(sapply(ENB, function(analysis) sapply(analysis, function(x) x$run_time)),
                    dim = c(n_analyses, n_sims),
                    dimnames = list("Analysis" = 0:(n_analyses - 1), "Simulation" = 1:n_sims))

  return(list(DEC = DEC, ENBP = ENB_PERFECT, ENBS = ENB_SAMPLE, RUN_TIME = RUN_TIME))
}
