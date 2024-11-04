# ValueAdapt

Install using *devtools::install_github("michaeldymock25/ValueAdapt")*.

This package contains the functions (with documentation and examples) to implement the *value-driven adaptive design* including:

- evpi() to estimate the expected value of perfect information using the Monte Carlo appoximation method
- evppi() to estimate the expected value of partial perfect information using either the Monte Carlo or non-parametric approximation methods
- evsi() to estimate the expected value of sample information using either the Monte Carlo, non-parametric or moment matching approximation methods

The *value-driven adaptive design* is implemented for the respiratory syncytial virus (RSV) [case study](./Examples/RSV). Note that the package will need to be loaded in order to run the case study and associated simulations (e.g., to load the evsi() function).