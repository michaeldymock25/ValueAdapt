
#' @title sim_betabinomial_trial
#' @import parallel
#' @description Simulates a series of trials using the value-driven adaptive design for a beta-binomial model
#' @param D Number of decision options
#' @param U Utility function that depends on the decision option and parameters
#' @param Theta Vector of true parameter values to use for simulation. Must be of length D.
#' @param t Named matrix of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H") for each analysis including the trial start-up decision. Columns are times and rows are analyses.
#' @param prop Vector containing current proportions of intervention use. Must sum to one.
#' @param cost Vector of costs to proceed between analyses. Must be of length n_analyses + 1 to include the trial start-up decision.
#' @param prior_par Matrix of parameters to specify the beta prior distribution with each row representing a different parameter contained in Theta.
#' @param n_analyses Maximum number of trial analyses where ENBS decisions are made (must be at least one)
#' @param n_samp Number of observations to sample between analyses for each decision option
#' @param n_sims Number of simulated trials to generate (must be at least two)
#' @param n_draws Number of parameter samples to draw from prior/posterior distributions
#' @param method Approximation method. Either MC for Monte-Carlo or NP for non-parametric (default). The moment matching method is not currently supported.
#' @param J Number of inner Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param K Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.
#' @param n_cores Number of cores to use if running in parallel
#' @return Matrix containing the expected value of sample information computed at each analysis for each simulated trial
#' @examples
#' D <- 2
#' U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
#' Theta <- c(0.5, 0.6)
#' n_analyses <- 2
#' t <- matrix(c(0, 1, 15,
#'               1, 2, 15,
#'               2, 3, 15),
#'             nrow = n_analyses + 1, ncol = 3, byrow = TRUE,
#'             dimnames = list(NULL, c("C", "A", "H")))
#' prop <- rep(1/D, D)
#' cost <- rep(0, n_analyses + 1)
#' prior_par <- matrix(rep(c(2, 3), D), nrow = D, byrow = TRUE)
#' sim_betabinomial_trial(D, U, Theta, t, prop, cost, prior_par,
#'                        n_analyses = n_analyses, n_samp = 50,
#'                        n_sims = 4)
#' @rdname sim_betabinomial_trial
#' @export
sim_betabinomial_trial <- function(D, U, Theta, t, prop, cost, prior_par, n_analyses = 1, n_samp = 1, n_sims = 1,
                                   n_draws = 10000, method = "NP", J = 10000, K = 10000, n_cores = 1){

  if(length(Theta) != D) stop("Theta must be a vector of length D")
  if(!(method %in% c("MC", "NP"))) stop("Method must be specified as MC or NP. Note that the moment matching method is not currently supported")
  if(nrow(t) != n_analyses + 1) stop("Matrix t must have rows for the number of analyses plus one")
  if(length(cost) != n_analyses + 1) stop("Cost vector must be of length the number of analyses plus one")

  dat_trial <- rbinom(n_sims*n_analyses*D, size = n_samp, Theta)
  dat_trial <- array(dat_trial,
                     dim = c(D, n_analyses, n_sims),
                     dimnames = list("Decision" = paste0("d", 1:D), "Analysis" = 1:n_analyses, "Simulation" = 1:n_sims))
  dat_trial <- apply(dat_trial, c("Decision", "Simulation"), cumsum)

  prior <- mclapply(1:n_sims,
                    function(sim){
                      lapply(1:D, function(d) rbeta(n_draws, prior_par[d,1], prior_par[d,2]))
                    },
                    mc.cores = n_cores)
  prior <- array(unlist(prior),
                 dim = c(n_draws, D, n_sims),
                 dimnames = list("Draw" = 1:n_draws, "Decision" = paste0("d", 1:D), "Simulation" = 1:n_sims))

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
                     dimnames = list("Draw" = 1:n_draws, "Decision" = paste0("d", 1:D), "Analysis" = 1:n_analyses, "Simulation" = 1:n_sims))

  post_fun <- function(J, x){
    matrix(unlist(lapply(1:D, function(d) rbeta(J, prior_par[d,1] + x[d], prior_par[d,2] + n_samp - x[d]))),
           nrow = J, ncol = D, dimnames = list(NULL, paste0("d", 1:D)))
  }
  samp_fun <- function(Theta){
    matrix(unlist(lapply(1:D, function(d) rbinom(n_draws, size = n_samp, prob = Theta[,d]))),
           nrow = n_draws, ncol = D, dimnames = list(NULL, paste0("d", 1:D)))
  }

  ENB_SAMPLE <- mclapply(1:n_sims,
                         function(sim){
                           lapply(0:n_analyses, function(analysis){
                              if(analysis == 0){
                                enb_sample(D = D, U = U, Theta = prior[,,sim], t = t[analysis+1,],
                                           prop = prop, cost = cost[analysis+1], method = method, post_fun = post_fun,
                                           samp_fun = samp_fun, J = J, K = K, stat_fun = function(x) x,
                                           model = paste0("s(", paste0("d", 1:D), ")", collapse = " + "))
                              } else {
                                enb_sample(D = D, U = U, Theta = posterior[,,analysis,sim], t = t[analysis+1,],
                                           prop = prop, cost = cost[analysis+1], method = method, post_fun = post_fun,
                                           samp_fun = samp_fun, J = J, K = K, stat_fun = function(x) x,
                                           model = paste0("s(", paste0("d", 1:D), ")", collapse = " + "))
                              }
                           })
                         },
                         mc.cores = n_cores)

  ENB_SAMPLE <- array(unlist(ENB_SAMPLE),
                      dim = c(n_analyses + 1, n_sims),
                      dimnames = list("Analysis" = 0:n_analyses, "Simulation" = 1:n_sims))

  return(ENB_SAMPLE)
}