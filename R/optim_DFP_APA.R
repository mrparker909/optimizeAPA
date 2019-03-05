# goal: implement DFP algorithm with arbitrary precision arithmetic
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
#' @title optim_DFP_APA
#' @description optim_DFP_APA arbitrary precision implementation of DFP (Davidon-Fletcher-Powel) quasi-newton optimization algorithm. The function func must take precBits as an arguement, to inform Rmpfr of the precision.
#' @param starts starting values for function parameters
#' @param func   function to optimize, first argument is a vector containing the parameters to be optimized over. Must also take as argument precBits.
#' @param precBits precision bits for determining stopping condition (how close to zero is considered small enough for the gradient). Note that larger precBits may require much larger maxSteps and lineSearchMaxSteps.
#' @param maxSteps maximum number of iterations for the optimization algorithm
#' @param lineSearchMaxSteps maximum number of iterations for each line search (occurs for every iteration of the optimization algorithm).
#' @param keepValues if TRUE will return all visited values during the optimization, rather than only the final values.
#' @param VERBOSE integer number, if greater than 0 will provide some console output during the optimization process
#' @param ... extra parameters passed on to func
#' @export
#' @examples 
#' # (1D example, compared against stats::optim)
#' optim_DFP_APA(15, func=function(x, precBits=64){(Rmpfr::mpfr(x,precBits)-10)^2}, maxSteps = 1000, precBits = 64)
#' optim(par = 15,     fn = function(x){(x-10)^2}, hessian = TRUE,   method="BFGS")
#' 
#' # (1D example, compared against stats::optim)
#' fun <- function(par, xdat, precBits=64) {
#'   l <- -1
#'   for(x in xdat) {
#'     l <- l * dpois_APA(x = x, lambda = par, precBits) 
#'   }
#'   return(l)
#' }
#' 
#' fun2 <- function(par, xdat) {
#'  -1*prod(dpois(x = xdat, lambda = par))
#' }
#' 
#' xdat <- c(8,11)#rpois(n = 2, lambda = 10)
#' starts <- 10
#' 
#' optim_DFP_APA(starts, fun, xdat=xdat, precBits=64)
#' optim(par = starts, fn = fun2, hessian = TRUE, method="BFGS", xdat=xdat)
#' 
#' # (2D example, compared against stats::optim)
#' fun2D <- function(par, xdat, ydat, precBits=64) {
#'   par <- exp(Rmpfr::mpfr(par, precBits))
#'   -1*(sum(log(dpois_APA(x = xdat, lambda = par[1], precBits)))+sum(log(dpois_APA(x = ydat, lambda = par[2], precBits))))
#' }
#' fun2D2 <- function(par, xdat, ydat) {
#'   par <- exp(par)
#'   -1*(sum(log(dpois(x = xdat, lambda = par[1])))+sum(log(dpois(x = ydat, lambda = par[2]))))
#' }

#' xdat2D <- c(1,2,3)
#' ydat2D <- c(5,8,9)
#' starts2D <- log(c(5,7))
#' 
#' trueValues <- c(log(mean(xdat2D)), log(mean(ydat2D)))
#' op1 <- optim_DFP_APA(starts2D, fun2D, xdat=xdat2D, ydat=ydat2D, precBits=64)
#' op2 <- optim(par = starts2D, fn = fun2D2, hessian = TRUE, method="BFGS", xdat=xdat2D, ydat=ydat2D)
optim_DFP_APA <- function(starts, func, precBits = 64, maxSteps=100, lineSearchMaxSteps=500, keepValues=FALSE, VERBOSE=0, ...) {
  require(Rmpfr)
  # 0) initial guess
  xk  <- Rmpfr::mpfr(starts, precBits=precBits)
  Ix  <- Rmpfr::mpfrArray(diag(rep(1, times=length(xk))), precBits=precBits, dim = c(length(xk),length(xk)))
  Bk  <- Ix
  iBk <- Bk
  
  # initialize lists
  x_list  <- updateList(new_el = xk)
  f_list  <- updateList(new_el = func(xk, precBits=precBits, ...))
  g_list  <- updateList(new_el = grad_FD_APA(func = func, x_val = xk, precBits=precBits, VERBOSE=VERBOSE, ...))
  
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
    ls <- lineSearch_APA(x_curr = x_list[[1]], 
                          dk = p_list[[1]], 
                          func = func,
                          grad_Fx = g_list[[1]],
                          precBits = precBits,
                          stepMod = steps,
                          lineSearchMaxSteps = lineSearchMaxSteps, ...)
    f_next <- ls$f_next
    x_next <- ls$x_next
    
    if(VERBOSE >= 2) print(paste("f_next =", format(f_next)))
    if(VERBOSE >= 2) print(paste("x_next =", format(x_next)))
    
    # 3) upadate lists
    f_list <- updateList(f_next, f_list, M = 1+length(f_list))
    if(VERBOSE >= 2) if(is.nan(f_list[[1]])) warning("WARNING: f_list[[1]] is NaN")
    x_list <- updateList(x_next, x_list, M = 1+length(x_list))
    if(VERBOSE >= 2) if(is.nan(t(x_list[[1]])%*%x_list[[1]])) warning("WARNING: x_list[[1]] is NaN")
    if(VERBOSE >= 2) print(paste("x_curr =", format(x_list[[1]])))
    s_list <- updateList(x_list[[1]]-x_list[[2]], s_list, M = 1+length(s_list))
    if(VERBOSE >= 2) if(is.nan(t(s_list[[1]])%*%s_list[[1]])) warning("WARNING: s_list[[1]] is NaN")
    
    # calculate gradient:
    g_next <- grad_FD_APA(func = func, x_val = x_next, precBits=precBits, stepMod=steps, VERBOSE=VERBOSE, ...)
    g_list <- updateList(g_next, g_list, M = 1+length(g_list))
    if(VERBOSE >= 2) if(is.nan(t(g_list[[1]])%*%g_list[[1]])) warning("WARNING: g_list[[1]] is NaN")
    y_list <- updateList(g_list[[1]]-g_list[[2]], y_list, M = 1+length(y_list))
    if(VERBOSE >= 2) if(is.nan(t(y_list[[1]])%*%y_list[[1]])) warning("WARNING: y_list[[1]] is NaN")
    
    if(sqrt(t(y_list[[1]])%*%y_list[[1]]) < Rmpfr::mpfr(0.5, precBits)^(precBits-1)) {
      warning("WARNING: difference in gradients is smaller than precision. Maybe use higher precision by increasing precBits?")
      converged=TRUE
    }
    # 4) calculate iBk, inverse approximate Hessian
    # DFP:
    y <- y_list[[1]]
    s <- s_list[[1]]
    iB <- iB_list[[1]]
    
    iB_next <- NULL
    if(t(y)%*%y==0 | t(s)%*%s==0) {
      iB_next <- iB_list[[1]]
    } else {
      iB_next <- iB + (s %*% t(s)) / c(t(s) %*% y) - (iB %*% y %*% t(y) %*% iB) / c(t(y) %*% iB %*% y)
    }
    
    if(any(is.na(iB_next))) { 
      # iB_next <- iB_list[[1]]
      print(paste("g = ", format(g_list[[1]])))
      print(paste("y = ", format(y)))
      print(paste("s = ", format(s)))
      print(paste("iB = ", format(iB)))
      stop(paste("WARNING: inverse hessian had NaN at iteration ", steps))
    }
    # update inverse Hessian list
    iB_list <- updateList(new_el = iB_next, iB_list, M = 1+length(iB_list))
    
    # check for convergence
    if( sqrt(g_list[[1]] %*% g_list[[1]]) < Rmpfr::mpfr(0.5,precBits)^(precBits-1)) {
      if(VERBOSE>=1) {
        print("CONVERGED:")
        print(paste("  stopping condition: ", format(sqrt(g_list[[1]] %*% g_list[[1]]))))
        print(paste("  was less than tolerance: ", format(Rmpfr::mpfr(0.5,precBits)^(precBits-4))))
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
    warning("WARNING: exceeded maxSteps, is maxSteps too small?")
    if(keepValues) {
      return(list(x=x_list, f=f_list, grad=g_list, inv_Hessian=iB_list, steps=steps, converged=converged))
    } else {
      return(list(x=x_list[[1]], f=f_list[[1]], grad=g_list[[1]], inv_Hessian=iB_list[[1]], steps=steps, converged=converged))
    }
  }
  
  if(keepValues) {
    return(list(x=x_list, f=f_list, grad=g_list, inv_Hessian=iB_list, steps=steps, converged=converged))
  } else {
    return(list(x=x_list[[1]], f=f_list[[1]], grad=g_list[[1]], inv_Hessian=iB_list[[1]], steps=steps, converged=converged))
  }
}


