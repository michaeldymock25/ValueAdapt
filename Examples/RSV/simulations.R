
## this script contains the simulations required to run the RSV case study (RSV.Rmd)
## the functions from the ValueAdapt package must be loaded

library(ValueAdapt)
library(here)

set.seed(12856421)

D <- 2
U <- function(d, Theta, t_1, t_2, lambda = 5200){
  cost <- c(590, 330)
  constant <- 1/1000000*sum(1.05^(1-(t_1:t_2)))*300000
  constant*(lambda*(1 - Theta[,d]) - cost[d])
}
t <- matrix(c(0, 2, 15,
              2, 3, 15,
              3, 4, 15),
            nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(NULL, c("C", "A", "H")))
prop <- c(1,0)
cost <- c(1 + 0.002*1000, 0.002*500, 0.002*500)
prior_par <- matrix(rep(c(4, 10), D), nrow = D, byrow = TRUE)
n_analyses <- 2
n_samp <- 250
n_sims <- 10000
n_cores <- 10

sims <- list(`Scenario 1` = sim_betabinomial_trial(D = D, U = U, Theta = c(0.10, 0.18), t = t, prop = prop, cost = cost,
                                                   prior_par = prior_par, n_analyses = n_analyses, n_samp = n_samp,
                                                   n_sims = n_sims, n_cores = n_cores),
             `Scenario 2` = sim_betabinomial_trial(D = D, U = U, Theta = c(0.10, 0.12), t = t, prop = prop, cost = cost,
                                                   prior_par = prior_par, n_analyses = n_analyses, n_samp = n_samp,
                                                   n_sims = n_sims, n_cores = n_cores))

sims <- data.frame(ENBS = unlist(sims),
                   Scenario = rep(paste("Scenario", 1:2), each = (n_analyses + 1)*n_sims),
                   Analysis = rep(paste("Analysis", 0:n_analyses), 2*n_sims),
                   Simulation = rep(rep(paste("Simulation", 1:n_sims), each = n_analyses + 1), 2))
rownames(sims) <- NULL

saveRDS(sims, here("Examples/RSV/simulations.rds"))