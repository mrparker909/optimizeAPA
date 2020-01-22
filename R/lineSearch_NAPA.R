#' @title lineSearch_NAPA
#' @description line searching algorithm used inside non-APA optimization algorithms.
#' @param x_curr starting x value
#' @param func function to perform line search on
#' @param dk direction vector to search along
#' @param grad_Fx (defaults to NULL) if not NULL, the gradient of func at x_curr
#' @param lineSearchMaxSteps maximum number of iterations before stopping
#' @param ... extra parameters passed on to func
#' @examples 
#' # simple 1D quadratic function optimization
#' lineSearch_NAPA(x_curr = 1, dk = -0.2, lineSearchMaxSteps = 100, func = function(x) {x^2})
#'
#' # simple 2D quadratic function optimization
#' func2 <- function(par, centerx=0, centery=0) {
#'   (par[1]-centerx)^2+(par[2]-centery)^2
#' }
#'
#' lineSearch_NAPA(x_curr = c(.015,-.015), dk = c(-0.2,0.2), func=func2)
#' 
#' lineSearch_NAPA(x_curr = c(2.05), dk = c(-0.2), func=function(x) { -prod(dpois(c(1,2,3), x))})
#' 
lineSearch_NAPA <- function(x_curr, dk, func, grad_Fx=NULL, stepMod=0, tolerance=10^-8, lineSearchMaxSteps = 10, ...) {
  
  delta <- 0.5#*2^-log(1+stepMod)
  alpha <- 0.05
  if(is.nan(t(dk)%*%dk)) warning("WARNING: direction vector dk is NaN")
  
  if(t(dk)%*%dk < 1) {
    dk <- dk / sqrt(sum(dk^2))
  }
  
  t <- 1
  f_curr <- func(x_curr, ...)
  x_best <- x_curr
  f_best <- f_curr
  
  x_next <- x_curr + dk
  f_next <- func(x_next, ...)
  if(f_next < f_best) {
    f_best <- f_next
    x_best <- x_next
  }
  
  gg <- NULL
  if(is.null(grad_Fx)) gg <- abs(grad_FD_NAPA(func=func, x_val=x_curr, stepMod, ...) %*% dk)
  else gg <- abs(grad_Fx %*% dk)
  
  # if(sqrt(gg) < tolerance) {
  #   return(list(x_next=x_best, f_next=f_best, iterations=0))
  # }
  # consider adding
  # if(f_next < f_curr) {
  #   return(list(x_next=x_next, f_next=f_next, iterations=0))
  # }
  ##
  lineSearchSteps <- 0
  lineSearching = TRUE
  if(f_next < f_curr - gg/2) {
    lineSearching = FALSE
  } else {
    while(lineSearching) {
      lineSearchSteps <- lineSearchSteps + 1
      if(lineSearchSteps > lineSearchMaxSteps) {
        # if(gg > tolerance) {
        #   warning("WARNING: exceeded lineSearchMaxSteps, is lineSearchMaxSteps too small?")
        # }
        if(f_next > f_curr) stop("ERROR: lineSearch returned larger function value")
        return(list(x_next=x_best, f_next=f_best, iterations=lineSearchSteps))
      }
      if(is.nan(f_next)) warning("WARNING: f_next is NaN")
      if(is.nan(f_curr)) warning("WARNING: f_curr is NaN")
      if(is.nan(gg)) warning("WARNING: gg is NaN")
      if(f_next < f_curr - alpha * t * gg) {
        lineSearching = FALSE
      } else {
        mod <- max(10^-4,2^-log(lineSearchSteps+stepMod))
        t <- delta*t*mod
        x_next <- x_curr + t * dk
        f_next <- func(x_next, ...)
        if(f_next < f_best) {
          f_best <- f_next
          x_best <- x_next
        }
      }
    }
  }
  if(f_next > f_curr) stop("ERROR: lineSearch returned larger function value")
  return(list(x_next=x_best, f_next=f_best, iterations=lineSearchSteps))
}




