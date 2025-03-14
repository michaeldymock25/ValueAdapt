% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ENB_PARTIAL_PERFECT.R
\name{enb_partial_perfect}
\alias{enb_partial_perfect}
\title{enb_partial_perfect}
\usage{
enb_partial_perfect(
  D,
  U,
  Theta_int,
  Theta_rem,
  t,
  prop,
  cost,
  method = "NP",
  K = 10000,
  cond_args = list(),
  cond_fun = NULL,
  model = NULL
)
}
\arguments{
\item{D}{Number of decision options}

\item{U}{Utility function that depends on the decision option, parameters and decision times t_1 and t_2}

\item{Theta_int}{Named matrix of parameter draws from prior/posterior distribution for parameters of interest}

\item{Theta_rem}{Named matrix of parameter draws from prior/posterior distribution for the remaining parameters (not of interest)}

\item{t}{Named vector of ascending decision times in years including the current time ("C"), analysis time ("A") and the time horizon ("H")}

\item{prop}{Vector containing current proportions of intervention use. Must sum to one.}

\item{cost}{Cost of sampling}

\item{method}{Approximation method. Either MC for Monte-Carlo, NP for non-parametric (default) or MM for moment matching. The moment matching method requires the evppi function to be run in advance using the non-parametric method to generate INB_partial.}

\item{K}{Number of outer Monte Carlo loops. Only required for the Monte Carlo approximation method.}

\item{cond_args}{List of arguments to be passed to cond_fun(). Defaults to an empty list.}

\item{cond_fun}{A function that generates sets of Theta_rem conditional on the current values of Theta_int. The only argument should be a list containing the parameter draws and additional parameters. Only required for the Monte Carlo approximation method.}

\item{model}{Generalised additive regression model specification (formula). Only required for the non-parametric approximation method. Should be a function of the parameters of interest.}
}
\value{
Expected net benefit of collecting partial perfect information
}
\description{
Computes the expected net benefit of collecting partial perfect information using either a Monte-Carlo or non-parametric approximation method
}
\examples{
# two parameters (one parameter of interest), two decision options
D <- 2
U <- function(d, Theta_int, Theta_rem, t_1, t_2)
  sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int + (d == 2)*Theta_rem)
N <- 10000
Theta_int <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_A"))
Theta_rem <- matrix(rbeta(N, 2, 3), nrow = N, ncol = 1, dimnames = list(NULL, "theta_B"))
t <- c(C = 0, A = 1, H = 15)
prop <- rep(1/D, D)
cost <- 0
cond_fun <- function(args) rbeta(args$N, 2, 3)
enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "MC",
                    cond_fun = cond_fun)
enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "NP",
                    model = "s(theta_A)")

# three parameters (two parameters of interest), three decision options
D <- 3
U <- function(d, Theta_int, Theta_rem, t_1, t_2)
  sum(1.05^(1-(t_1:t_2)))*((d == 1)*Theta_int[,"theta_A"] +
                           (d == 2)*Theta_int[,"theta_B"] +
                           (d == 3)*Theta_rem)
N <- 10000
Theta <- matrix(c(rbeta(N, 2, 3), rbeta(N, 2, 3), rbeta(N, 2, 3)),
                nrow = N, ncol = 3, dimnames = list(NULL, c("theta_A", "theta_B", "theta_C")))
Theta_int <- Theta[,c("theta_A", "theta_B")]
Theta_rem <- matrix(Theta[,"theta_C"], nrow = N, ncol = 1, dimnames = list(NULL, "theta_C"))
t <- c(C = 0, A = 1, H = 15)
prop <- rep(1/D, D)
cost <- 0
cond_fun <- function(args) rbeta(args$N, 2, 3)
enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "MC",
                    cond_fun = cond_fun)
enb_partial_perfect(D, U, Theta_int, Theta_rem, t, prop, cost, method = "NP",
                    model = "te(theta_A, theta_B)")
}
