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
lineSearch_APA <- function(x_curr, dk, func, grad_Fx=NULL, precBits=64, stepMod=0, lineSearchMaxSteps = 100, ...) {
  ten <- Rmpfr::mpfr(10, precBits)
  two <- Rmpfr::mpfr(2, precBits)
  delta <- Rmpfr::mpfr(0.5, precBits)#*2^-log(1+stepMod)#*ten^(-stepMod)
  alpha <- Rmpfr::mpfr(0.05, precBits)#*2^-log(1+stepMod)#*ten^(-1-(stepMod))
  x_curr <- Rmpfr::mpfr(x_curr, precBits) 
  
  dk <- Rmpfr::mpfr(dk, precBits)
  if(t(dk)%*%dk < 1) {
    dk <- Rmpfr::mpfr(dk, precBits)/sqrt(sum(dk^2))
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
    gg <- abs(grad_FD_APA(func  = func,
                          x_val = x_curr,
                          precBits = precBits, ...) %*% dk)
  } else {
    grad_Fx <- Rmpfr::mpfr(grad_Fx, precBits)
    gg <- abs(grad_Fx %*% dk)
  }
  
  # if(gg < two^-precBits) {
  #   return(list(x_next=x_best, f_next=f_best, iterations=0))
  # }
  
  # consider adding this?
  # if(f_next < f_curr) {
  #   if(abs(f_next-f_curr) > gg/2) return(list(x_next=x_next, f_next=f_next, iterations=0))
  # }
  #if(f_next < f_curr) { return(list(x_next=x_next, f_next=f_next, iterations=0)) }
  ##
  
  lineSearchSteps <- 0
  lineSearching = TRUE
  while(lineSearching) {
    lineSearchSteps <- lineSearchSteps + 1
    if(lineSearchSteps > lineSearchMaxSteps) {
      if(gg > two^(-precBits/2)) {
        warning("WARNING: exceeded lineSearchMaxSteps, is lineSearchMaxSteps too small?")
      }
      if(f_next > f_curr) stop("ERROR: lineSearch returned larger function value")
      return(list(x_next=x_best, f_next=f_best, iterations=lineSearchSteps))
    }
    if(any(is.na(f_next))) stop("ERROR: f_next is NaN")
    if(any(is.na(f_curr))) stop("ERROR: f_curr is NaN")
    if(any(is.na(gg))) stop("ERROR: gg is NaN")
    if(f_next <  f_curr - alpha * t * gg) {
      lineSearching = FALSE
    } else {
      #t <- delta*2^(-stepMod)*t
      mod <- max(two^-(precBits/8),two^-log(lineSearchSteps+stepMod))
      t <- delta*t*2^-log(lineSearchSteps) #delta*ten^(-lineSearchSteps+1)*t
      x_next <- x_curr + t * dk
      f_next <- func(x_next, precBits=precBits, ...)
      if(f_next < f_best) {
        f_best <- f_next
        x_best <- x_next
      }
    }
  }
  return(list(x_next=x_best, f_next=f_best, iterations=lineSearchSteps))
}




