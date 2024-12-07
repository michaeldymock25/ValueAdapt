% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ENB_SAMPLE.R
\name{enb_sample}
\alias{enb_sample}
\title{enb_sample}
\usage{
enb_sample(
  D,
  U,
  Theta,
  t,
  prop,
  cost,
  method = "NP",
  J = 10000,
  K = 10000,
  samp_fun = NULL,
  post_fun = NULL,
  stat_fun = NULL,
  model = NULL,
  INB_partial = NULL,
  Q = 50
)
}
\arguments{
\item{D}{Number of decision options}

\item{U}{Utility function that depends on the decision option, parameters and decision times t_1 and t_2}

\item{Theta}{Named matrix of parameter draws from prior/posterior distribution}

\item{t}{Named vector of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H")}

\item{prop}{Vector containing current proportions of intervention use. Must sum to one.}

\item{cost}{Cost of sampling}

\item{method}{Approximation method. Either MC for Monte-Carlo, NP for non-parametric (default) or MM for moment matching. The moment matching method requires the evppi function to be run in advance using the non-parametric method to generate INB_partial.}

\item{J}{Number of inner Monte Carlo loops. Only required for the Monte Carlo approximation method.}

\item{K}{Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.}

\item{samp_fun}{A function that generates a sample of data based on the parameter draws (posterior predictive draw). The only argument is the matrix of parameter draws Theta. Returns a matrix of data.}

\item{post_fun}{A function that generates draws from the posterior distribution given the previous distribution and a sample of data generated by the samp_fun function. The first argument must be the number of draws to generate and the second argument must be the sample of data generated by the samp_fun function. Only required for the Monte Carlo and moment matching approximation methods.}

\item{stat_fun}{A function to generate a low-dimensional summary statistic based on the sample of data generated by the samp_fun function. The only argument is a vector of sampled data. Only required for the non-parametric approximation method.}

\item{model}{Generalised additive regression model specification (formula). Only required for the non-parametric approximation method.}

\item{INB_partial}{Samples of INB for the parameters of interest generated from the evppi function using the non-parametric approximation method. Only required for the moment matching approximation method.}

\item{Q}{Number of model reruns to estimate the expected variance of the posterior net benefit. Only required for the moment matching approximation method.}
}
\value{
Expected net benefit of collecting sample information
}
\description{
Computes the expected net benefit of collecting sample information using either a Monte-Carlo, non-parametric or moment matching approximation method
}
\examples{
# one parameter, two decision options
D <- 2
U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*(-1)^(d-1)*(Theta - 0.4)
N <- 10000
Theta <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta"))
t <- c(C = 0, A = 1, H = 15)
prop <- rep(1/D, D)
cost <- 0
samp_fun <- function(Theta) matrix(rbinom(nrow(Theta), size = 10, prob = Theta),
                                   nrow = nrow(Theta), ncol = 1, dimnames = list(NULL, "theta"))
post_fun <- function(J, x) rbeta(J, 2 + x, 3 + 10 - x)
stat_fun <- function(x) x/10
enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
           model = "s(theta)")

# two parameters (one parameter of interest), two decision options
D <- 2
U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
N <- 10000
Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3)),
                nrow = N, ncol = 2, dimnames = list(NULL, c("theta_A", "theta_B")))
t <- c(C = 0, A = 1, H = 15)
prop <- rep(1/D, D)
cost <- 0
samp_fun <- function(Theta) matrix(rbinom(nrow(Theta), size = 10, prob = Theta[,"theta_A"]),
                                   nrow = nrow(Theta), ncol = 1, dimnames = list(NULL, "theta_A"))
post_fun <- function(J, x) matrix(c(rbeta(J, 2 + x, 3 + 10 - x), rbeta(J, 2, 3)),
                                  nrow = J, ncol = 2,
                                  dimnames = list(NULL, c("theta_A", "theta_B")))
stat_fun <- function(x) x/10
enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
           model = "s(theta_A)")
U_enbppi <- function(d, Theta_int, Theta_rem, t_1, t_2)
  sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int + (d == 2)*Theta_rem)
Theta_int <- matrix(Theta[,"theta_A"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
Theta_rem <- matrix(Theta[,"theta_B"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
out <- enb_partial_perfect(D, U_enbppi, Theta_int, Theta_rem, t, prop, cost, method = "NP",
                           model = "s(theta_A)")
enb_sample(D, U, Theta, t, prop, cost, method = "MM", samp_fun = samp_fun, post_fun = post_fun,
           INB_partial = out$INB_partial)

# three parameters (two parameters of interest), three decision options
D <- 3
U <- function(d, Theta, t_1, t_2) sum(1.05^(1-(t_1:t_2)))*Theta[,d]
N <- 10000
Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
                nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
t <- c(C = 0, A = 1, H = 15)
prop <- rep(1/D, D)
cost <- 0
samp_fun <- function(Theta) apply(Theta[,c("theta_A", "theta_B")], 2,
                                  function(theta) rbinom(length(theta), size = 10, prob = theta))
post_fun <- function(J, x) matrix(c(rbeta(J, 2 + x["theta_A"], 3 + 10 - x["theta_A"]),
                                    rbeta(J, 2 + x["theta_B"], 3 + 10 - x["theta_B"]),
                                    rbeta(J, 2, 3)),
                                  nrow = J, ncol = 3,
                                  dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
stat_fun <- function(x) x/10
enb_sample(D, U, Theta, t, prop, cost, method = "MC", samp_fun = samp_fun, post_fun = post_fun)
enb_sample(D, U, Theta, t, prop, cost, method = "NP", samp_fun = samp_fun, stat_fun = stat_fun,
           model = "te(theta_A, theta_B)")
U_enbppi <- function(d, Theta_int, Theta_rem, t_1, t_2)
  sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int[,"theta_A"] +
                             (d == 2)*Theta_int[,"theta_B"] +
                             (d == 3)*Theta_rem)
Theta_int <- Theta[,c("theta_A", "theta_B")]
Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
cond_fun <- function(J, Theta_int) rbeta(J, 2, 3)
out <- enb_partial_perfect(D, U_enbppi, Theta_int, Theta_rem, t, prop, cost, method = "NP",
                           cond_fun = cond_fun, model = "te(theta_A, theta_B)")
enb_sample(D, U, Theta, t, prop, cost, method = "MM", samp_fun = samp_fun,
           post_fun = post_fun, INB_partial = out$INB_partial)
}
