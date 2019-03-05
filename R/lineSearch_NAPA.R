#' @description line searching algorithm used inside non-APA optimization algorithms.
#' @param x_curr starting x value
#' @param func function to perform line search on
#' @param dk direction vector to search along
#' @param grad_Fx (defaults to NULL) if not NULL, the gradient of func at x_curr
#' @param lineSearchMaxSteps (defaults to 100) maximum number of iterations before stopping
#' @param ... extra parameters passed on to func
#' @examples 
#' # simple 1D quadratic function optimization
#' lineSearch_NAPA(x_curr = 1, dk = -0.2, lineSearchMaxSteps = 100, func = function(x) {x^2})
#'
#' # simple 2D quadratic function optimization
#' func2 <- function(par, centerx=0, centery=0) {
#'   (par[1]-centerx)^2+(par[2]-centery)^2
#' }
# '
#' lineSearch_NAPA(x_curr = c(.015,-.015), dk = c(-0.2,0.2), func=func2)
lineSearch_NAPA <- function(x_curr, dk, func, grad_Fx=NULL, lineSearchMaxSteps = 100, ...) {
  delta <- 0.5
  if(is.nan(t(dk)%*%dk)) warning("WARNING: direction vector dk is NaN")
  
  t <- 1
  f_curr <- func(x_curr, ...)
  x_next <- x_curr + dk
  f_next <- func(x_next, ...)
  gg <- NULL
  if(is.null(grad_Fx)) gg <- abs(grad_FD(func=func, x_val=x_curr, ...) %*% dk)
  else gg <- grad_Fx %*% dk
    
  lineSearchSteps <- 0
  lineSearching = TRUE
  while(lineSearching) {
    lineSearchSteps <- lineSearchSteps + 1
    if(lineSearchSteps > lineSearchMaxSteps) {
      warning("WARNING: exceeded lineSearchMaxSteps, is lineSearchMaxSteps too small?")
      return(list(x_next=x_next, f_next=f_next, iterations=lineSearchSteps))
    }
    if(is.nan(f_next)) warning("WARNING: f_next is NaN")
    if(is.nan(f_curr)) warning("WARNING: f_curr is NaN")
    if(is.nan(gg)) warning("WARNING: gg is NaN")
    if(f_next < f_curr - 0.5 * t * gg) {
      lineSearching = FALSE
    } else {
      t <- delta*t
      x_next <- x_curr + t * dk
      f_next <- func(x_next, ...)
    }
  }
  return(list(x_next=x_next, f_next=f_next, iterations=lineSearchSteps))
}




