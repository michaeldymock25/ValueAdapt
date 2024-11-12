# ValueAdapt

Install using *devtools::install_github("michaeldymock25/ValueAdapt")*. If the *EVSI* package is not installed then it will need to be installed first using *devtools::install_github("annaheath/EVSI")*. This package requires an R version >= 3.5.0 if the examples code is to be used.

This package contains the functions (with documentation and examples) to implement the *value-driven adaptive design* including:

- evpi() to estimate the expected value of perfect information using a Monte Carlo appoximation method
- enb_perfect() to estimate the expected net benefit of collecting perfect information using a Monte Carlo appoximation method
- evppi() to estimate the expected net benefit of collecting partial perfect information using either a Monte Carlo or non-parametric approximation method
- enb_partial_perfect() to estimate the expected value of partial perfect information for a trial using either a Monte Carlo or non-parametric approximation method
- evsi() to estimate the expected value of sample information using either a Monte Carlo, non-parametric or moment matching approximation method
- enb_sample() to estimate the expected net benefit of collecting sample information using either a Monte Carlo, non-parametric or moment matching approximation method
- sim_betabinomial_trial() to simulate a series of trials using the value-driven adaptive design for a beta-binomial model

The enb_*() functions extend the base functions by incorporating the length of the trial into the value of information computation via an ascending vector containing the current time, time of analysis and time horizon. The time of analysis influence the value of the information collected (e.g., information that takes longer to collect will be less influential on the value of information). The functions return the expected net benefit of collecting information (including the cost of data collection).

The *value-driven adaptive design* is implemented for the respiratory syncytial virus (RSV) [case study](./Examples/RSV). Note that the package will need to be loaded in order to run the case study and associated simulations (e.g., to load the enb_sample() function).