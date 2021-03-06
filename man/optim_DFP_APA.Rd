% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim_DFP_APA.R
\name{optim_DFP_APA}
\alias{optim_DFP_APA}
\title{optim_DFP_APA}
\usage{
optim_DFP_APA(
  starts,
  func,
  tolerance = 10^-10,
  precBits = 64,
  maxSteps = 100,
  lineSearchMaxSteps = 100,
  keepValues = FALSE,
  Memory = 100,
  outFile = NULL,
  ...
)
}
\arguments{
\item{starts}{vector of starting values for the function parameters}

\item{func}{function to optimize, first argument is a vector containing the parameters to be optimized over. Must also take as argument precBits.}

\item{tolerance}{tolerance for determining stopping condition (how close to zero is considered small enough for the gradient). Note that smaller tolerance may require much larger maxSteps and lineSearchMaxSteps. Tolerance should either be larger than 10^-10, or an arbitrary precision number, such as: tolerance=Rmpfr::mpfr(10^-20, precBits=128)}

\item{precBits}{determines number of bits of precision for all numbers and calculations, for standard double precision use precBits=53.}

\item{maxSteps}{maximum number of iterations for the optimization algorithm}

\item{lineSearchMaxSteps}{maximum number of iterations for each line search (occurs for every iteration of the optimization algorithm).}

\item{keepValues}{if TRUE will return all visited values (up to the number specified by Memory) during the optimization, rather than only the final values.}

\item{Memory}{maximum number of iterations to remember (default is 100)}

\item{outFile}{if not NULL, name of file to save results to (will be overwritten at each iteration).}

\item{...}{extra parameters passed on to func}
}
\description{
optim_DFP_APA is an arbitrary precision implementation of DFP (Davidon-Fletcher-Powel) quasi-newton optimization algorithm. The function to be optimized, func, must take precBits as an arguement, to inform Rmpfr of the precision being used.
}
\examples{
# (1D example, compared against stats::optim)
optim_DFP_APA(15, func=function(x, precBits=64){(Rmpfr::mpfr(x,precBits)-10)^2}, maxSteps = 1000, precBits = 64)
optim(par = 15, fn = function(x){(x-10)^2}, hessian = TRUE, method="BFGS")

# (1D example, compared against stats::optim)
fun <- function(par, xdat, precBits=64) {
  l <- -1
  for(x in xdat) {
    l <- l * dpois_APA(x = x, lambda = par, precBits) 
  }
  return(l)
}

fun2 <- function(par, xdat) {
 -1*prod(dpois(x = xdat, lambda = par))
}

xdat <- c(8,11)#rpois(n = 2, lambda = 10)
starts <- 10

optim_DFP_APA(starts, fun, xdat=xdat, precBits=64)
optim(par = starts, fn = fun2, hessian = TRUE, method="BFGS", xdat=xdat)

# (2D example, compared against stats::optim)
fun2D <- function(par, xdat, ydat, precBits=64) {
  par <- exp(Rmpfr::mpfr(par, precBits))
  -1*(sum(log(dpois_APA(x = xdat, lambda = par[1], precBits)))+sum(log(dpois_APA(x = ydat, lambda = par[2], precBits))))
}
fun2D2 <- function(par, xdat, ydat) {
  par <- exp(par)
  -1*(sum(log(dpois(x = xdat, lambda = par[1])))+sum(log(dpois(x = ydat, lambda = par[2]))))
}
xdat2D <- c(1,2,3)
ydat2D <- c(5,8,9)
starts2D <- log(c(5,7))

trueValues <- c(log(mean(xdat2D)), log(mean(ydat2D)))
op1 <- optim_DFP_APA(starts2D, fun2D, xdat=xdat2D, ydat=ydat2D, precBits=64, keepValues=T)
op2 <- optim(par = starts2D, fn = fun2D2, hessian = TRUE, method="BFGS", xdat=xdat2D, ydat=ydat2D)

plotConvergence(op1, digits=20)
  
# N dimensional quadratic
funcND <- function(par, center, precBits=64) {
  par <- Rmpfr::mpfr(par, precBits)
  sum((par-center)^2)
} 
funcND2 <- function(par, center) {
  sum((par-center)^2)
} 

op1 <- optim_DFP_APA(starts = c(0,0,0,0,0), func = funcND, center=c(1,2,3,4,5), keepValues=T, Memory=10)
optim(par = c(0,0,0,0,0), fn = funcND2, hessian = TRUE, method="BFGS", center=c(1,2,3,4,5))

plotConvergence(op1)


funcexpND <- function(par, center, precBits=64) {
  par <- Rmpfr::mpfr(par, precBits)
  log(sum(exp((par-center)^2)))
} 

op1 <- optim_DFP_APA(starts = c(0,0,0,0,0), func = funcexpND, center=c(1,-2,3,-4,5), keepValues=T, Memory=100)
plotConvergence(op1)
}
