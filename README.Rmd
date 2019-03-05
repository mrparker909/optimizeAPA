# optimizeAPA
Arbitrary Precision Arithmetic (APA) optimization algorithms implemented in R.

Throughout I will use the abreviations APA = "arbitrary precision arithmetic", and NAPA="non arbitrary precision arithmetic".

The APA in this package is achieved using the packages Rmpfr and gmp.

# Installation

```{r}
install.packages("devtools")
devtools::install_github("mrparker909/optimizeAPA")
```

## Currently Implemented:

- optim_DFP_NAPA (non-APA version of DFP algorithm)
- optim_DFP_APA (APA version of DFP algorithm)

## Currently In Progress:

- optim_BFGS_NAPA (non-APA version of BFGS algorithm)
- optim_BFGS_APA (APA version of BFGS algorithm)

# Some Illustrative Examples

## optimDFP_APA

```{r}
library(optimizeAPA)

# notice that the function must be made using APA arithmetic from the package Rmpfr
# the function inputs are first: the parameter to be optimized over (x)
# and must also include the variable precBits, to be passed to all APA numbers.
# Any other function parameters (such as center) can be passed to the function
# when calling optim_DFP_APA

# function we would like to optimize:
quadratic <- function(x, center, precBits=64) {
  x_apa <- Rmpfr::mpfr(x,precBits) # convert input to APA
  (x_apa-center)^2
}

# APA optimization
optim_DFP_APA(starts=0, func=quadratic, center=10, precBits = 64)
```
The minimum is attained at $x=10$, where the function is found to be $f(x_min)=0$, and the gradient is $grad(f(x_min)) = 0$.

We can compare these results to the optim_DFP_NAPA results, and to the regular optim results: 

```{r}
# function we would like to optimize:
quadratic_NAPA <- function(x, center) {
  (x-center)^2
}

# NAPA optimization
optim_DFP_NAPA(starts=0, func=quadratic_NAPA, center=10)
optim(par=0, fn=quadratic_NAPA, hessian=TRUE, method="BFGS", center=10)
```

It can be desirable to retain the convergence path for diagnosis and visualization purposes. This can be done easily with both the APA and NAPA version of the optimization algorithm:



