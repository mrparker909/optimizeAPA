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

#' @title carryForwardNA
#' @description Checks if front of list is NA, then replaces it with the most recent value which is not NA.
#' @param thelist the list to check for NA head
#' @export
carryForwardNA <- function(thelist) {
  if(all(is.na(unlist(thelist)))) { stop("ERROR: every element of list is NA, no carry forward possible") }

  l = length(thelist)
  i = 2
  while(any(is.na(thelist)) & i <= l) {
    if(!any(is.na(thelist[[i]])) &
       any(is.na(thelist[[i-1]]))) {
      thelist[[i-1]] <- thelist[[i]]
    }
    i = i+1
    if(i>l) { i=2 }
  }
  if(any(is.na(thelist[[1]]))) { stop("ERROR: NA in list position 1, maybe increase precBits") }
  return(thelist)
}

# calculate gradient using finite difference method (requires 2*N function evaluations, where N is the dimension of x)
grad_FD_NAPA <- function(func, x_val, stepMod=0, stepsize = .Machine$double.eps^(1/3), ...) {
  stepsize <- stepsize/(1+stepMod)
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
grad_FD_APA <- function(func, x_val, stepMod=0, precBits=64, ...) {
  if(class(x_val)!="mpfr") 
    { x <- Rmpfr::mpfr(x_val, precBits=precBits) } 
  else 
    { x <- x_val }
  
  stepsize <- Rmpfr::mpfr(2,precBits=precBits)^(-precBits/2-stepMod+2)
  len_x <- length(x)
  delta_f <- list()
  for(i in seq(len_x)) {
    step <- Rmpfr::mpfr(rep(0, times = len_x), precBits=precBits)
    step[i] <- stepsize
    delta_f[[i]] <- func(x+step, precBits=precBits, ...) - func(x-step, precBits=precBits, ...)
  }
  
  delta_f <- methods::new("mpfr", unlist(delta_f))/(2*stepsize)
  return(delta_f)
}

#' @title Arbitrary Precision Poisson density function.
#' @description Poisson density function with arbitrary precision. Uses library Rmpfr for arbitrary precision arithmetic.
#' @param x quantile
#' @param lambda poisson mean parameter
#' @param precBits number of bits of precision
#' @export
#' @examples 
#' # compare apa to non apa dpois:
#' sum(dpois_APA(x = 0:70, lambda = 10, precBits = 128))
#' sum(dpois(x = 0:40, lambda = 10)) # sums to 1, indicating that dpois(x=41, lambda=10) should be zero
#' dpois(x=41, lambda=10) # is not zero
#' dpois_APA(x=41, lambda=10, precBits = 128)
dpois_APA <- function(x, lambda, precBits=128) {
  ind <- which(x < 0)
  x[ind] <- 0
  ind2 <- which(!gmp::is.whole(x))
  x[ind2] <- 0
  # create arbitrary precision arithmetic numbers
  if(class(x)!="mpfr") { x_apa <- Rmpfr::mpfr(x, precBits = precBits) }
  else { x_apa <- x }
  if(class(lambda)!="mpfr") { lambda_apa <- Rmpfr::mpfr(lambda, precBits = precBits) }
  else{ lambda_apa <- lambda}
  # calculate the density
  dens <- exp(-1*lambda_apa)*lambda_apa^(x_apa)/Rmpfr::mpfr(gmp::factorialZ(x),precBits)
  dens[ind]  <- 0
  dens[ind2] <- 0
  
  return(dens)
}




#' @title Arbitrary Precision Logistic density function.
#' @description Logistic density function with arbitrary precision. Uses library Rmpfr for arbitrary precision arithmetic.
#' @param x quantile
#' @param location Location parameter (default 0)
#' @param scale    Scale parameter (default 1)
#' @param precBits Number of bits of precision
#' @export
plogis_APA <- function(x, location = 0, scale=1, precBits=128) {
  # create APA numbers
  if(class(x)!="mpfr") { x_apa <- (location-Rmpfr::mpfr(x, precBits=precBits))/scale }
  else { x_apa <- (location-x)/scale }
  
  # calculate density
  # F(x) = 1 / (1 + exp(-(x-m)/s))
  dens <- 1 / (1+exp(x_apa))
  return(dens)
}


#' @title Arbitrary Precision Binomial density function.
#' @description Binomial density function with arbitrary precision. Uses library Rmpfr for arbitrary precision arithmetic.
#' @param x quantile (number of Bernoulli successes with success probability prob)
#' @param size number of Bernoulli trials to perform with success probability prob
#' @param prob probability of success for each independent Bernoulli trial
#' @param precBits number of bits of precision
#' @export
#' @examples 
#' # compare apa to non apa dpois:
#' sum(dbinom_APA(x = 0:10, size = 10, prob=0.2, precBits = 256))
#' sum(dbinom(x = 0:10, size = 10, prob=0.2))
dbinom_APA <- function(x, size, prob, precBits=128) {
  if(class(x)!="mpfr") { x_apa <- Rmpfr::mpfr(x, precBits=precBits) }
  else{ x_apa <- x }
  if(class(size)!="mpfr") { size_apa <- Rmpfr::mpfr(size, precBits=precBits) }
  else { size_apa <- size }
  if(class(prob)!="mpfr") { prob_apa <- Rmpfr::mpfr(prob, precBits=precBits) }
  else{ prob_apa <- prob }
  q_apa <- Rmpfr::mpfr(1, precBits=precBits)-prob_apa
  
  dens <- gmp::chooseZ(size,x) * (prob_apa)^(x_apa) * (q_apa)^(size-x_apa)
  dens <- Rmpfr::mpfr(dens, precBits=precBits)
  return(dens)
}




#' @title plotConvergence
#' @description Plot the path to convergence of an optimization algorithm (such as optim_DFP_APA).
#' @param optimOutput The result of running an optimization algorithm with keepValues=TRUE
#' @param variable Default is NULL, will print convergence plot for all variables, otherwise an integer vector (such as c(1,3)) indicating which variables to print convergence plot for.
#' @param digits Default is 5. Number of digits to print on convergence plot.
#' @param labels If TRUE, print function values as labels on the plot.
#' @examples 
#' # N dimensional quadratic
#' funcND <- function(par, center, precBits=64) {
#'   par <- Rmpfr::mpfr(par, precBits)
#'   log(sum(exp((par-center)^2)))
#' } 
#' opt <- optimizeAPA::optim_DFP_APA(starts = c(0,0,0), func = funcND, center=c(1,2,3), keepValues=TRUE)
#' plotConvergence(opt, variable=c(1,3), labels=F)
plotConvergence <- function(optimOutput, variable=NULL, digits=5, labels=T) {
  require(ggplot2)
  op <- optimOutput
  len <- length(op$f)
  lenx <- length(op$grad[[1]])
  
  dat <- data.frame()
  dat_arrows <- data.frame()
  fv <- sapply(X = op$f, FUN = as.numeric)
  for(i in 1:len) {
    for(j in 1:lenx) {
      if(!is.null(variable)) {
        if(!(j %in% variable)) {
          next()
        } 
      }
      xvij <- as.numeric(op$x[[i]][[j]])
      dat <- rbind(dat,data.frame(par=j, xv=xvij, FunctionValue=fv[i], steps=(op$steps-i+1)))
    }
  }
  
  xv_final <- unlist(lapply(op$x[[1]], FUN = as.numeric ))
  if(!is.null(variable)) {
    xv_final <- xv_final[variable] 
  }

  
  dat_text <- data.frame(
    label = paste("Final parameter value:",format(xv_final,digits=digits)),
    par   = unique(dat$par)
  )
  
  p <- ggplot(data=dat, aes(group=par)) + 
    geom_path(aes(x=xv, y=steps),
              arrow = arrow(length=unit(0.03,"npc"), angle=40, ends="first", type = "closed"), color="forestgreen", size=1.25) +
    geom_point(aes(x=xv, y=steps, color=FunctionValue), size=2.75) +
    facet_wrap(.~par, scales = "free") +
    geom_text(data = dat_text,
      mapping = aes(x = -Inf, y = op$steps+1, label = label),
      hjust   = -0.5) +
    xlab("Visited x Values") + ylab("Algorithm Step") +
    ggtitle("Convergence Path") + theme_classic() + scale_y_continuous(breaks=1:op$steps)
  
  if(labels) {
    p <- p + geom_label(aes(label=format(FunctionValue,digits=digits), x=xv, y=steps+0.5), size=2, alpha=0.6)
  }
  
  plot(p)
}

