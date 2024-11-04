
## this script contains the simulations required to run the RSV case study (RSV.Rmd)

library(parallel)
library(ValueAdapt)
library(here)

cores <- 10
seed <- 356723556

## simulate trial data

N_scenario <- 2
N_arm <- 2
N_sim <- 10000
N_analysis <- 4
N_prop <- 250
dat_trial <- list(`Scenario 1` = array(rbinom(N_sim*N_analysis*N_arm, size = N_prop, c(0.10, 0.18)),
                                       dim = c(N_arm, N_analysis, N_sim),
                                       dimnames = list(c("II", "MV"),
                                                       paste("Analysis", 1:N_analysis),
                                                       paste("Simulation", 1:N_sim))),
                  `Scenario 2` = array(rbinom(N_sim*N_analysis*N_arm, size = N_prop, c(0.10, 0.10)),
                                       dim = c(N_arm, N_analysis, N_sim),
                                       dimnames = list(c("II", "MV"),
                                                       paste("Analysis", 1:N_analysis),
                                                       paste("Simulation", 1:N_sim))))
dat_trial <- lapply(dat_trial, function(dat) apply(dat, c(1,3), cumsum))

## generate posterior distributions

N_draw <- 10000
posterior <- lapply(dat_trial, function(dat)
                mclapply(1:N_sim,
                         function(sim){
                           lapply(1:N_analysis, function(analysis)
                              lapply(1:N_arm, function(arm)
                                 rbeta(N_draw, 4 + dat[analysis, arm, sim], 20 + N_prop*analysis - dat[analysis, arm, sim])))
                         },
                         mc.cores = cores,
                         mc.set.seed = seed))
posterior <- list(`Scenario 1` = array(unlist(posterior$`Scenario 1`),
                                       dim = c(N_draw, N_arm, N_analysis, N_sim),
                                       dimnames = list(paste("Draw", 1:N_draw),
                                                       c("II", "MV"),
                                                       paste("Analysis", 1:N_analysis),
                                                       paste("Simulation", 1:N_sim))),
                  `Scenario 2` = array(unlist(posterior$`Scenario 2`),
                                       dim = c(N_draw, N_arm, N_analysis, N_sim),
                                       dimnames = list(paste("Draw", 1:N_draw),
                                                       c("II", "MV"),
                                                       paste("Analysis", 1:N_analysis),
                                                       paste("Simulation", 1:N_sim))))

## define utility function, sampling function, statistic function and model

U <- function(d, Theta, lambda = 5200){
  cost <- ifelse(d == 1, 590, 330)
  1/1000000*sum(1.05^(1-(1:15)))*300000*(-lambda*Theta[,d] - cost)
}

samp_fun <- function(Theta, N_prop = 250)
  matrix(c(rbinom(nrow(Theta), size = N_prop, prob = Theta[,"II"]),
           rbinom(nrow(Theta), size = N_prop, prob = Theta[,"MV"])),
         nrow = nrow(Theta), ncol = 2, dimnames = list(NULL, c("II", "MV")))

stat_fun <- function(x) x

model <- "s(II) + s(MV)"

EVSI <- lapply(posterior, function(dat)
           mclapply(1:N_sim,
                    function(sim){
                      lapply(1:N_analysis, function(analysis)
                        evsi(D = 2, U = U, Theta = dat[,,analysis,sim], samp_fun = samp_fun, stat_fun = stat_fun, model = model))
                    },
                    mc.cores = cores,
                    mc.set.seed = seed))
EVSI <- array(unlist(EVSI),
              dim = c(N_analysis, N_sim, N_scenario),
              dimnames = list(paste("Analysis", 1:N_analysis),
                              paste("Simulation", 1:N_sim),
                              paste("Scenario", 1:N_scenario)))

sims <- data.frame(EVSI = as.vector(EVSI),
                   Scenario = rep(paste("Scenario", 1:N_scenario), each = N_sim*N_analysis),
                   Simulation = rep(rep(paste("Simulation", 1:N_sim), each = N_analysis), N_scenario),
                   Analysis = rep(paste("Analysis", 1:N_analysis), N_scenario*N_sim))

saveRDS(sims, here("Examples/RSV/simulations.rds"))