#' @title lineSearch_APA
#' @description line searching algorithm used inside APA optimization algorithms.
#' @param x_curr starting x value
#' @param dk direction vector to search along
#' @param func function to perform line search on
#' @param grad_Fx (defaults to NULL) if not NULL, the gradient of func at x_curr
#' @param precBits bits of precision
#' @param lineSearchMaxSteps (defaults to 100) maximum number of iterations before stopping
#' @param ... extra parameters passed on to func
#' @examples 
#' # simple 1D quadratic function optimization
#' lineSearch_APA(x_curr = 1, dk = -0.2, lineSearchMaxSteps = 100, func = function(x, precBits) {Rmpfr::mpfr(x,precBits)^2})
#' lineSearch_APA(x_curr = 1, dk = -0.2, lineSearchMaxSteps = 100, func = function(x, precBits) {Rmpfr::mpfr(x-1,precBits)^2})
#'
#' # simple 2D quadratic function optimization
#' func2 <- function(par, centerx=0, centery=0, precBits=64) {
#'   par <- Rmpfr::mpfr(par, precBits)
#'   (par[1]-centerx)^2+(par[2]-centery)^2
#' }
#'
#' lineSearch_APA(x_curr = c(.015,-.015), dk = c(-0.2,0.2), func=func2)
#' 
#' lineSearch_APA(x_curr = c(2.06), dk = c(-0.2), func=function(x, precBits) { 
#'  l <- -1
#'  for(xi in c(1,2,3)) {
#'    l <- l * dpois_APA(xi, x, precBits)
#'  }
#'  return(l)
#'  })
lineSearch_APA <- function(x_curr, dk, func, grad_Fx=NULL, precBits=64, stepMod=0, lineSearchMaxSteps = 100, ...) {
  
  two   <- Rmpfr::mpfr(2, precBits)
  delta <- Rmpfr::mpfr(0.5, precBits)
  alpha <- Rmpfr::mpfr(0.05, precBits)
  if(class(x_curr)!="mpfr") { x_curr <- Rmpfr::mpfr(x_curr, precBits) }
  if(class(dk)!="mpfr") { dk <- Rmpfr::mpfr(dk, precBits) }
  
  if(t(dk)%*%dk < 1) {
    dk <- dk/sqrt(sum(dk^2))
  }
  
  t <- 1
  f_curr <- func(x_curr, precBits=precBits, ...)
  x_best <- x_curr
  f_best <- f_curr
  
  x_next <- x_curr + dk
  f_next <- func(x_next, precBits=precBits, ...)
  if(f_next < f_best) {
    f_best <- f_next
    x_best <- x_next
  }
  
  gg <- NULL
  if(is.null(grad_Fx)) {
    gg <- abs(grad_FD_APA(func     = func,
                          x_val    = x_curr,
                          precBits = precBits, 
                          stepMod  = stepMod, ...) %*% dk)
  } else {
    if(class(grad_Fx)!="mpfr") { grad_Fx <- Rmpfr::mpfr(grad_Fx, precBits) }
    gg <- abs(grad_Fx %*% dk)
  }
  
  
  lineSearchSteps <- 0
  lineSearching = TRUE
  if(f_next < f_curr - gg/2) {
    lineSearching = FALSE
  } else {
    while(lineSearching) {
      lineSearchSteps <- lineSearchSteps + 1
      if(lineSearchSteps > lineSearchMaxSteps) {
        # if(gg > two^(-precBits/2)) {
        #   warning("WARNING: exceeded lineSearchMaxSteps, is lineSearchMaxSteps too small?")
        # }
        if(f_next > f_curr) stop("ERROR: lineSearch returned larger function value")
        return(list(x_next=x_best, f_next=f_best, iterations=lineSearchSteps))
      }
      if(any(is.na(f_next))) stop("ERROR: f_next is NA")
      if(any(is.na(f_curr))) stop("ERROR: f_curr is NA")
      if(any(is.na(gg))) stop("ERROR: gg is NA")
      if(f_next <  f_curr - alpha * t * gg) {
        lineSearching = FALSE
      } else {
        mod <- two^-log(lineSearchSteps+stepMod)
        t <- delta*t*mod
        x_next <- x_curr + t * dk
        f_next <- func(x_next, precBits=precBits, ...)
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




