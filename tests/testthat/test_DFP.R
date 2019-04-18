context("Testing DFP Algorithms")

test_that("optim_DFP_NAPA", {
  op1 <- optim_DFP_NAPA(15, func=function(x){(x-10)^2}, tolerance = 10^-16)
  expect_equal(as.numeric(abs(op1$x-10)) < 10^-8, TRUE)
  
  op1 <- optim_DFP_NAPA(-15, func=function(x){(x-10)^2}, tolerance = 10^-16)
  op2 <- optim(par = -15, fn = function(x){(x-10)^2}, hessian = TRUE,   method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-10, TRUE)
  
  fun <- function(par, xdat) {
    -1*prod(dpois(x = xdat, lambda = par))
  }

  op1 <- optim_DFP_NAPA(10, fun, xdat=c(8,11), tolerance=10^-8)
  expect_equal(as.numeric(abs(op1$x-9.5)) < 10^-8, TRUE)
  
  op1 <- optim_DFP_NAPA(8, fun, xdat=c(8,11), tolerance=10^-10)
  expect_equal(as.numeric(abs(op1$x-9.5)) < 10^-8, TRUE)
  
  fun2D <- function(par, xdat, ydat) {
    par <- exp(par)
    -1*(sum(dpois(x = xdat, lambda = par[1], log=TRUE))+sum(dpois(x = ydat, lambda = par[2], log=TRUE)))
  }
  
  xdat2D <- c(1,2,3)
  ydat2D <- c(5,8,9)
  starts2D <- log(c(5,7))
   
  op1 <- optim_DFP_NAPA(starts = starts2D, func = fun2D, xdat=xdat2D, ydat=ydat2D, tolerance=10^-6)
  op2 <- optim(par = starts2D, fn = fun2D, method="BFGS", xdat=xdat2D, ydat=ydat2D)
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-4, c(TRUE,TRUE))
  
  
  funND <- function(par, centers) {
    sum((par-centers)^6)
  }

  par     <- c(1,2,3,4,5)
  centers <- c(5,4,3,2,1)

  op1 <- optim_DFP_NAPA(starts = par, func = funND, centers=centers, tolerance = 10^-8)
  expect_equal(as.numeric(abs(op1$x-centers)) < 10^-1, c(TRUE,TRUE,TRUE,TRUE,TRUE))
})



test_that("optim_DFP_APA", {
  op1 <- optim_DFP_APA(15, func=function(x, precBits=53){(Rmpfr::mpfr(x,precBits)-10)^2})
  op2 <- optim(par = 15, fn = function(x){(x-10)^2}, method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-10, TRUE)
  
  op1 <- optim_DFP_APA(-15, func=function(x, precBits=53){(Rmpfr::mpfr(x,precBits)-10)^2})
  op2 <- optim(par = -15, fn = function(x){(x-10)^2}, method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-10, TRUE)

  fun <- function(par, xdat, precBits=53) {
    l <- -1
    for(x in xdat) {
      l <- l * dpois_APA(x = x, lambda = par, precBits) 
    }
    return(l)
  }

  op1 <- optim_DFP_APA(starts = 10, func = fun, xdat=c(8,11), precBits=128, lineSearchMaxSteps = 200)
  expect_equal(as.numeric(abs(op1$x-mean(c(8,11)))) < 10^-6, TRUE)

  
  funND <- function(par, centers, precBits=53) {
    par     <- Rmpfr::mpfr(par, precBits)
    centers <- Rmpfr::mpfr(centers, precBits)

    sum((par-centers)^6)
  }

  precBits <- 128
  par     <- c(1,2,3,4,5)
  centers <- c(5,4,3,2,1)

  op1 <- optim_DFP_APA(starts = par, func = funND, centers=centers, precBits=precBits, maxSteps = 150, lineSearchMaxSteps = 500)
  expect_equal(as.numeric(abs(op1$x-centers)) < 10^-1, c(TRUE,TRUE,TRUE,TRUE,TRUE))
})
