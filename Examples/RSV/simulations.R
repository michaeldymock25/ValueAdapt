
## this script contains the simulations required to run the RSV case study (RSV.Rmd)

library(ValueAdapt)
library(here)

set.seed(356723556)

D <- 2
U <- function(d, Theta, lambda = 5200){
  cost <- ifelse(d == 1, 590, 330)
  1/1000000*sum(1.05^(1-(1:15)))*300000*(-lambda*Theta[,d] - cost)
}
prior <- matrix(rep(c(4, 10), D), nrow = D, byrow = TRUE)
n_analyses <- 4
n_samp <- 250
n_sims <- 10000
n_cores <- 10

sims <- list(`Scenario 1` = sim_betabinomial_trial(D = D, U = U, Theta = c(0.10, 0.18), prior = prior, n_analyses = n_analyses,
                                                   n_samp = n_samp, n_sims = n_sims, n_cores = n_cores),
             `Scenario 2` = sim_betabinomial_trial(D = D, U = U, Theta = c(0.10, 0.10), prior = prior, n_analyses = n_analyses,
                                                   n_samp = n_samp, n_sims = n_sims, n_cores = n_cores))

sims <- data.frame(EVSI = unlist(sims),
                   Scenario = rep(paste("Scenario", 1:2), each = n_analyses*n_sims),
                   Analysis = rep(paste("Analysis", 1:n_analyses), 2*n_sims),
                   Simulation = rep(rep(paste("Simulation", 1:n_sims), each = n_analyses), 2))
rownames(sims) <- NULL

saveRDS(sims, here("Examples/RSV/simulations.rds"))