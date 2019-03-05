#' @title lineSearch_APA
#' @description line searching algorithm used inside APA optimization algorithms.
#' @param x_curr starting x value
#' @param func function to perform line search on
#' @param dk direction vector to search along
#' @param grad_Fx (defaults to NULL) if not NULL, the gradient of func at x_curr
#' @param lineSearchMaxSteps (defaults to 100) maximum number of iterations before stopping
#' @param precBits bits of precision
#' @param ... extra parameters passed on to func
#' @examples 
#' # simple 1D quadratic function optimization
#' lineSearch_APA(x_curr = 1, dk = -0.2, lineSearchMaxSteps = 100, func = function(x, precBits) {Rmpfr::mpfr(x,precBits)^2})
#'
#' # simple 2D quadratic function optimization
#' func2 <- function(par, centerx=0, centery=0, precBits=64) {
#'   par <- Rmpfr::mpfr(par, precBits)
#'   (par[1]-centerx)^2+(par[2]-centery)^2
#' }
# '
#' lineSearch_APA(x_curr = c(.015,-.015), dk = c(-0.2,0.2), func=func2)
#' 
#' lineSearch_APA(x_curr = c(2.06), dk = c(-0.2), func=function(x, precBits) { 
#'  l <- -1
#'  for(xi in c(1,2,3)) {
#'    l <- l * dpois_APA(xi, x, precBits)
#'  }
#'  return(l)
#'  })
lineSearch_APA <- function(x_curr, dk, func, grad_Fx=NULL, precBits=64, stepMod=0, lineSearchMaxSteps = 100, VERBOSE=0, ...) {
  require(Rmpfr)
  ten <- Rmpfr::mpfr(10, precBits)
  delta <- Rmpfr::mpfr(0.5, precBits)*ten^(-stepMod)
  alpha <- Rmpfr::mpfr(0.5, precBits)*ten^(-1-log(stepMod))
  x_curr <- Rmpfr::mpfr(x_curr, precBits) 
  dk <- Rmpfr::mpfr(dk, precBits)
  
  
  if(is.nan(t(dk)%*%dk)) warning("WARNING: direction vector dk is NaN")
  if(any(is.na(grad_Fx))) stop("grad_Fx is NaN")
  t <- 1
  f_curr <- func(x_curr, precBits=precBits, ...)
  x_next <- x_curr + dk
  f_next <- func(x_next, precBits=precBits, ...)
  gg <- NULL
  if(is.null(grad_Fx)) {gg <- abs(grad_FD_APA(func=func, x_val=x_curr, precBits=precBits, ...) %*% dk)}
  else {
    grad_Fx <- Rmpfr::mpfr(grad_Fx, precBits)
    gg <- grad_Fx %*% dk
  }
  if(VERBOSE >=3) print(format(gg))
  
  lineSearchSteps <- 0
  lineSearching = TRUE
  while(lineSearching) {
    lineSearchSteps <- lineSearchSteps + 1
    if(lineSearchSteps > lineSearchMaxSteps) {
      warning("WARNING: exceeded lineSearchMaxSteps, is lineSearchMaxSteps too small?")
      return(list(x_next=x_next, f_next=f_next, iterations=lineSearchSteps))
    }
    if(any(is.na(f_next))) stop("ERROR: f_next is NaN")
    if(any(is.na(f_curr))) stop("ERROR: f_curr is NaN")
    if(any(is.na(gg))) stop("ERROR: gg is NaN")
    if(f_next <  f_curr - alpha * t * gg | f_next - f_curr < Rmpfr::mpfr(0.5, precBits)^(precBits-1)) {
      lineSearching = FALSE
    } else {
      t <- delta*t
      x_next <- x_curr + t * dk
      f_next <- func(x_next, precBits=precBits, ...)
    }
  }
  return(list(x_next=x_next, f_next=f_next, iterations=lineSearchSteps))
}




