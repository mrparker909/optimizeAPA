#' @title updateList
#' @description Adds new element to top of list, removing bottom elements if new length is greater tham M
#' @param new_el   element to add to top of list
#' @param old_list list to append new_el to
#' @param M        maximum size of list (will truncate by removing bottom if length is greater than M)
updateList <- function(new_el, old_list = list(), M = 1) {
  old_list <- c(list(new_el), old_list) 
  if(length(old_list) > M) {
    old_list <- old_list[1:M]
  }
  return(old_list)
}


# calculate gradient using finite difference method (requires 2*N function evaluations, where N is the dimension of x)
grad_FD_NAPA <- function(func, x_val, stepsize = .Machine$double.eps^(1/3), ...) {
  #stepsize <- rep(stepsize, times=length(x))
  len_x <- length(x_val)
  delta_f <- numeric(len_x)
  for(i in seq(len_x)) {
    step <- rep(0, times = len_x)
    step[i] <- stepsize
    delta_f[i] <- func(x_val+step, ...) - func(x_val-step, ...)
  }
  
  return(delta_f/(2*stepsize))
}

# calculate gradient using finite difference method (requires 2*N function evaluations, where N is the dimension of x)
# grad_FD_APA(func = function(x, precBits){x^2}, x=0)
# grad_FD_APA(func = function(x, precBits){x^2}, x=1)
# grad_FD_APA(func = function(x, precBits){x^2}, x=2)
# grad_FD_APA(func = function(x, precBits){ x[1]^2 + x[2]^2}, x=c(1,3))
grad_FD_APA <- function(func, x_val, precBits=64, ...) {
  require(Rmpfr)
  x <- Rmpfr::mpfr(x_val, precBits=precBits)
  stepsize <- Rmpfr::mpfr(0.5, precBits=precBits)^(precBits-4) #Rmpfr::mpfr(.Machine$double.eps, precBits=precBits)^Rmpfr::mpfr(1/3), precBits=precBits)
  len_x <- length(x)
  delta_f <- list()
  for(i in seq(len_x)) {
    step <- Rmpfr::mpfr(rep(0, times = len_x), precBits=precBits)
    step[i] <- stepsize
    delta_f[[i]] <- func(x+step, precBits=precBits, ...) - func(x-step, precBits=precBits, ...)
  }
  
  delta_f <- new("mpfr", unlist(delta_f))/(Rmpfr::mpfr(2, precBits=precBits)*stepsize)
  return(delta_f)
}


#' @description Poisson density function with arbitrary precision. Uses library Rmpfr for arbitrary precision arithmetic.
#' @param x quantile
#' @param lambda poisson mean parameter
#' @param precBits number of bits of precision
#' @examples 
#' # compare apa to non apa dpois:
#' sum(dpois_apa(x = 0:100, lambda = 900, prec = 100))
#' sum(dpois(x = 0:100, lambda = 900))
dpois_apa <- function(x, lambda, precBits) {
  require(Rmpfr)
  # create arbitrary precision arithmetic numbers
  x_apa      <- Rmpfr::mpfr(x, precBits = precBits)
  lambda_apa <- Rmpfr::mpfr(lambda, precBits = precBits)
  one_apa    <- Rmpfr::mpfr(1,precBits=precBits)
  # calculate the density
  dens <- exp(-one_apa*lambda_apa)*lambda_apa^(x_apa)/factorial(x_apa)
  if(is.nan(dens)) {
    stop(paste0("DENSITY IS NaN:", "\n x_apa=",format(x_apa), "\n lambda_apa=",format(lambda_apa)))
  }
  return(dens)
}