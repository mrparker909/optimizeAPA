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