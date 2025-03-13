
## this script contains the correction computed for the RSV case study (RSV.Rmd)

library(ValueAdapt)
library(here)

set.seed(12856421)

num_cores <- 20

D <- 2
U <- function(d, Theta, t_1, t_2, lambda = 5200){
  cost <- c(590, 330)
  constant <- 1/1000000*sum(1.05^(1-(t_1:t_2)))*300000
  constant*(lambda*(1 - Theta[,d]) - cost[d])
}
prop <- c(1,0)
n_draws <- 10000
prior <- matrix(c(rbeta(n_draws, 4, 20), rbeta(n_draws, 4, 20)),
                nrow = n_draws, ncol = D,
                dimnames = list(NULL, c("p_II", "p_MV")))
dat_1 <- matrix(c(33, 47,
                  53, 87,
                  66, 125,
                  91, 166), nrow = 4, ncol = 2, byrow = TRUE, dimnames = list(paste("Analysis", 1:4), c("II", "MV")))
post_1 <- lapply(1:4, function(i)
  matrix(c(rbeta(n_draws, 4 + dat_1[i,"II"], 20 + 250*i - dat_1[i,"II"]),
           rbeta(n_draws, 4 + dat_1[i,"MV"], 20 + 250*i - dat_1[i,"MV"])),
         nrow = n_draws, ncol = 2,
         dimnames = list(NULL, c("p_II", "p_MV"))))
dat_2 <- matrix(c(27 , 26,
                  48 , 65,
                  82 , 99,
                  103, 121), nrow = 4, ncol = 2, byrow = TRUE, dimnames = list(paste("Analysis", 1:4), c("II", "MV")))
post_2 <- lapply(1:4, function(i)
  matrix(c(rbeta(n_draws, 4 + dat_2[i,"II"], 20 + 250*i - dat_2[i,"II"]),
           rbeta(n_draws, 4 + dat_2[i,"MV"], 20 + 250*i - dat_2[i,"MV"])),
         nrow = n_draws, ncol = 2,
         dimnames = list(NULL, c("p_II", "p_MV"))))
t <- matrix(c(0, 2, 15,
              2, 3, 15,
              3, 4, 15,
              4, 5, 15),
            nrow = 4, ncol = 3, byrow = TRUE,
            dimnames = list(NULL, c("C", "A", "H")))
delta <- 1
gamma <- 0.002
cost <- c(delta + gamma*1000, gamma*500, gamma*500, 0)
samp_fun <- function(args){
  c(II = rbinom(1, size = 250, prob = args[["p_II"]]),
    MV = rbinom(1, size = 250, prob = args[["p_MV"]]))
}
post_args <- list(n_prev = 0, II_prev = 0, MV_prev = 0)
post_args_1 <- list(list(n_prev = 250, II_prev = 33, MV_prev = 47),
                    list(n_prev = 500, II_prev = 53, MV_prev = 87))
post_args_2 <- list(list(n_prev = 250, II_prev = 27, MV_prev = 26),
                    list(n_prev = 500, II_prev = 48, MV_prev = 65))
post_fun <- function(args){
  matrix(c(rbeta(args[["N"]],
                 4 + (args[["II_prev"]] + args[["II"]]),
                 20 + (args[["n_prev"]] + 250) - (args[["II_prev"]] + args[["II"]])),
           rbeta(args[["N"]],
                 4 + (args[["MV_prev"]] + args[["MV"]]),
                 20 + (args[["n_prev"]] + 250) - (args[["MV_prev"]] + args[["MV"]]))),
         nrow = args[["N"]], ncol = 2, dimnames = list(NULL, c("p_II", "p_MV")))
}

## calculate ENBS_0

correct_0 <- list(n_analyses = 2,
                  t_update_args = list(t = t[1:2,]),
                  t_update = function(t_update_args) t_update_args[["t"]][t_update_args[["analysis"]],],
                  cost_update_args = list(cost = cost[1:2]),
                  cost_update = function(cost_update_args) cost_update_args[["cost"]][cost_update_args[["analysis"]]],
                  post_args_update = function(post_args_update_args)
                                       list(n_prev = post_args_update_args[["n_prev"]] + 250,
                                            II_prev = post_args_update_args[["II_prev"]] + post_args_update_args[["II"]],
                                            MV_prev = post_args_update_args[["MV_prev"]] + post_args_update_args[["MV"]]))

start_time_0 <- Sys.time()
enbs_0 <- enb_sample(D = D, U = U, Theta = prior, t = t[1,], prop = prop, cost = cost[1], method = "MC",
                     samp_fun = samp_fun, post_args = post_args, post_fun = post_fun, correct = correct_0, num_cores = num_cores)
end_time_0 <- Sys.time()
run_time_0 <- difftime(end_time_0, start_time_0, units = "mins")

print(enbs_0)
print(run_time_0)

## calculate ENBS_1 for each scenario

correct_1 <- list(n_analyses = 2,
                  t_update_args = list(t = t[2:3,]),
                  t_update = function(t_update_args) t_update_args[["t"]][t_update_args[["analysis"]],],
                  cost_update_args = list(cost = cost[2:3]),
                  cost_update = function(cost_update_args) cost_update_args[["cost"]][cost_update_args[["analysis"]]],
                  post_args_update = function(post_args_update_args)
                                       list(n_prev = post_args_update_args[["n_prev"]] + 250,
                                            II_prev = post_args_update_args[["II_prev"]] + post_args_update_args[["II"]],
                                            MV_prev = post_args_update_args[["MV_prev"]] + post_args_update_args[["MV"]]))

start_time_11 <- Sys.time()
enbs_11 <- enb_sample(D = D, U = U, Theta = post_1[[1]], t = t[2,], prop = prop, cost = cost[2], method = "MC",
                      samp_fun = samp_fun, post_args = post_args_1[[1]], post_fun = post_fun, correct = correct_1, num_cores = num_cores)
end_time_11 <- Sys.time()
run_time_11 <- difftime(end_time_11, start_time_11, units = "mins")

print(enbs_11)
print(run_time_11)

start_time_12 <- Sys.time()
enbs_12 <- enb_sample(D = D, U = U, Theta = post_2[[1]], t = t[2,], prop = prop, cost = cost[2], method = "MC",
                      samp_fun = samp_fun, post_args = post_args_2[[1]], post_fun = post_fun, correct = correct_1, num_cores = num_cores)
end_time_12 <- Sys.time()
run_time_12 <- difftime(end_time_12, start_time_12, units = "mins")

print(enbs_12)
print(run_time_12)

## calculate ENBS_2 for each scenario

correct_2 <- list(n_analyses = 2,
                  t_update_args = list(t = t[3:4,]),
                  t_update = function(t_update_args) t_update_args[["t"]][t_update_args[["analysis"]],],
                  cost_update_args = list(cost = cost[3:4]),
                  cost_update = function(cost_update_args) cost_update_args[["cost"]][cost_update_args[["analysis"]]],
                  post_args_update = function(post_args_update_args)
                                       list(n_prev = post_args_update_args[["n_prev"]] + 250,
                                            II_prev = post_args_update_args[["II_prev"]] + post_args_update_args[["II"]],
                                            MV_prev = post_args_update_args[["MV_prev"]] + post_args_update_args[["MV"]]))

start_time_21 <- Sys.time()
enbs_21 <- enb_sample(D = D, U = U, Theta = post_1[[2]], t = t[3,], prop = prop, cost = cost[3], method = "MC",
                      samp_fun = samp_fun, post_args = post_args_1[[2]], post_fun = post_fun, correct = correct_2, num_cores = num_cores)
end_time_21 <- Sys.time()
run_time_21 <- difftime(end_time_21, start_time_21, units = "mins")

print(enbs_21)
print(run_time_21)

start_time_22 <- Sys.time()
enbs_22 <- enb_sample(D = D, U = U, Theta = post_2[[2]], t = t[3,], prop = prop, cost = cost[3], method = "MC",
                      samp_fun = samp_fun, post_args = post_args_2[[2]], post_fun = post_fun, correct = correct_2, num_cores = num_cores)
end_time_22 <- Sys.time()
run_time_22 <- difftime(end_time_22, start_time_22, units = "mins")

print(enbs_22)
print(run_time_22)

out <- list(enbs_0      = enbs_0,
            run_time_0  = run_time_0,
            enbs_11     = enbs_11,
            run_time_11 = run_time_11,
            enbs_12     = enbs_12,
            run_time_12 = run_time_12,
            enbs_21     = enbs_21,
            run_time_21 = run_time_21,
            enbs_22     = enbs_22,
            run_time_22 = run_time_22)

saveRDS(out, here("Examples/RSV/RSV_correction.rds"))