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


# goal: implement DFP algorithm
# 0) initial guess xk = starts, initial hessian Bk = Ipxp
# 1) solve for pk: Bk * pk = -grad f(xk) (use inverse(Bk), Ipxp in first iteration, otherwise from step 4)
# 2) line search to find step size a ~= argmin_a f(xk + a*pk)
# 3) update values: set sk = ak * pk, update xnext = xk + sk
#    yk = grad f(xnext) - grad f(xk)
# 4) calculate inverse(Bnext) http://www.stat.cmu.edu/~ryantibs/convexopt-F13/lectures/11-QuasiNewton.pdf
#    p = delta x= x_next - x_curr, 
#    q = delta grad f(x) = grad f(x_next) - grad f(x_curr)
#    iB_next + p %*% t(p) / c(t(p) %*% q) - iB %*% q %*% t(q) %*% iB / c(t(q) %*% iB %*% q)
# 5) stop algorithm if ||grad (f(xnext))|| is close enough to zero, otherwise goto step 1
#' @description optim_DFP_NAPA non-arbitrary precision implementation of DFP (Davidon-Fletcher-Powel) quasi-newton optimization algorithm
#' @param starts starting values for function
#' @param func   function to optimize
#' @param tolerance tolerance for determining stopping condition (how close to zero is considered small enough for the gradient). Note that smaller tolerance may require much larger maxSteps and lineSearchMaxSteps.
#' @param maxSteps maximum number of iterations for the optimization algorithm
#' @param lineSearchMaxSteps maximum number of iterations for each line search (occurs for every iteration of the optimization algorithm).
#' @param keepValues if TRUE will return all visited values during the optimization, rather than only the final values.
#' @param VERBOSE integer number, if greater than 0 will provide some console output during the optimization process
#' @param ... extra parameters passed on to func
#' @export
#' @examples 
#' # (1D example, compared against stats::optim)
#' optim_DFP_NAPA(15, func=function(x){(x-10)^2}, maxSteps = 1000, tolerance = 10^-6)
#' optim(par = 15,     fn = function(x){(x-10)^2}, hessian = TRUE,   method="BFGS")
#' 
#' # (1D example, compared against stats::optim)
#' fun <- function(par, xdat) {
#'  -1*prod(dpois(x = xdat, lambda = par))
#' }
#' 
#' xdat <- c(8,11)#rpois(n = 2, lambda = 10)
#' starts <- 10
#' 
#' optim_DFP_NAPA(starts, fun, xdat=xdat, tolerance=10^-10)
#' optim(par = starts, fn = fun, hessian = TRUE, method="BFGS", xdat=xdat)
#' 
#' # (2D example, compared against stats::optim)
#' fun2D <- function(par, xdat, ydat) {
#'   par <- exp(par)
#'   -1*(sum(dpois(x = xdat, lambda = par[1], log=TRUE))+sum(dpois(x = ydat, lambda = par[2], log=TRUE)))
#' }
#'     
#' xdat2D <- c(1,2,3)
#' ydat2D <- c(5,8,9)
#' starts2D <- log(c(5,7))
#' 
#' op1 <- optim_DFP_NAPA(starts2D, fun2D, xdat=xdat2D, ydat=ydat2D, tolerance=10^-6)
#' op2 <- optim(par = starts2D, fn = fun2D, hessian = TRUE, method="BFGS", xdat=xdat2D, ydat=ydat2D)
optim_DFP_NAPA <- function(starts, func, tolerance = 10^-10, maxSteps=100, lineSearchMaxSteps=100, keepValues=FALSE, VERBOSE=0, ...) {
  # 0) initial guess
  xk  <- starts
  Ix  <- diag(rep(1, times=length(xk)))
  Bk  <- Ix
  iBk <- Bk
  
  # initialize lists
  x_list  <- updateList(new_el = xk)
  f_list  <- updateList(new_el = func(xk, ...))
  g_list  <- updateList(new_el = grad_FD_NAPA(func = func, x_val = xk, ...))
  
  # B_list  <- updateList(new_el = Bk)
  iB_list <- updateList(new_el = iBk)

  p_list <- list()
  s_list <- list()
  y_list <- list()
  
  # begin optimization
  converged <- FALSE
  steps     <- 0
  while(steps < maxSteps & !converged) {
    steps <- steps+1

    # 1) solve for pk (direction vector)
    pk <- -1 * iB_list[[1]] %*% g_list[[1]]
    p_list <- updateList(pk, p_list, M = 1+length(p_list))
    
    # 2) line search to find step size t ~= argmin_t f(xk + t*pk)
    ls <- lineSearch_NAPA(x_curr = x_list[[1]], 
                          dk = p_list[[1]], 
                          func = func,
                          grad_Fx = g_list[[1]],
                          lineSearchMaxSteps = lineSearchMaxSteps, ...)
    f_next <- ls$f_next
    x_next <- ls$x_next
    
    # 3) upadate lists
    f_list <- updateList(f_next, f_list, M = 1+length(f_list))
    if(VERBOSE >= 2) if(is.nan(f_list[[1]])) warning("WARNING: f_list[[1]] is NaN")
    x_list <- updateList(x_next, x_list, M = 1+length(x_list))
    if(VERBOSE >= 2) if(is.nan(t(x_list[[1]])%*%x_list[[1]])) warning("WARNING: x_list[[1]] is NaN")
    print(paste("x_curr=", x_list[[1]]))
    s_list <- updateList(x_list[[1]]-x_list[[2]], s_list, M = 1+length(s_list))
    if(VERBOSE >= 2) if(is.nan(t(s_list[[1]])%*%s_list[[1]])) warning("WARNING: s_list[[1]] is NaN")
    g_list <- updateList(grad_FD_NAPA(func = func, x_val = x_next, ...), g_list, M = 1+length(g_list))
    if(VERBOSE >= 2) if(is.nan(t(g_list[[1]])%*%g_list[[1]])) warning("WARNING: g_list[[1]] is NaN")
    y_list <- updateList(g_list[[1]]-g_list[[2]], y_list, M = 1+length(y_list))
    if(VERBOSE >= 2) if(is.nan(t(y_list[[1]])%*%y_list[[1]])) warning("WARNING: y_list[[1]] is NaN")
    
    if(sqrt(t(y_list[[1]])%*%y_list[[1]]) < .Machine$double.eps) warning("WARNING: difference in gradients is smaller than machine precision. Consider using optim_DFP_APA for higher precision.")
  
    # 4) calculate iBk, inverse approximate Hessian
    # DFP:
    y <- y_list[[1]]
    s <- s_list[[1]]
    iB <- iB_list[[1]]
    iB_next <- iB + s %*% t(s) / c(t(s) %*% y) - iB %*% y %*% t(y) %*% iB / c(t(y) %*% iB %*% y)
    
    if(any(is.na(iB_next))) { 
      iB_next <- iB_list[[1]]
      if(VERBOSE >= 1) warning(paste("WARNING: inverse hessian had NaN at iteration ", steps))
    }
    # update inverse Hessian list
    iB_list <- updateList(new_el = iB_next, iB_list, M = 1+length(iB_list))
    
    # check for convergence
    if( sqrt(g_list[[1]] %*% g_list[[1]]) < tolerance) {
      if(VERBOSE>=1) {
        print("CONVERGED:")
        print(paste("  stopping condition: ", sqrt(g_list[[1]] %*% g_list[[1]])))
        print(paste("  was less than tolerance: ", tolerance))
      }
      converged = TRUE
    }
    
    if(converged) {
      if(VERBOSE >= 2) {
        print(paste("lineSearchSteps = ", ls$iterations))
        print(paste("Steps = ", steps))
      }
    }
  }
  
  if(steps >= maxSteps) {
    stop("ERROR: exceeded maxSteps, is maxSteps too small?")
  }
  
  if(keepValues) {
    return(list(x=x_list, f=f_list, grad=g_list, inv_Hessian=iB_list, steps=steps, converged=converged))
  } else {
    return(list(x=x_list[[1]], f=f_list[[1]], grad=g_list[[1]], inv_Hessian=iB_list[[1]], steps=steps, converged=converged))
  }
}


