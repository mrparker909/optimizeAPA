context("Testing DFP Algorithms")

test_that("optim_DFP_NAPA", {
  op1 <- optim_DFP_NAPA(15, func=function(x){(x-10)^2}, tolerance = 10^-8)
  expect_equal(as.numeric(abs(op1$x-10)) < 10^-8, TRUE)
  
  op1 <- optim_DFP_NAPA(-15, func=function(x){(x-10)^2}, tolerance = 10^-8)
  op2 <- optim(par = -15, fn = function(x){(x-10)^2}, hessian = TRUE,   method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-8, TRUE)
  
  fun <- function(par, xdat) {
    -1*prod(dpois(x = xdat, lambda = par))
  }

  op1 <- optim_DFP_NAPA(10, fun, xdat=c(8,11), tolerance=10^-8, maxSteps = 300)
  expect_equal(as.numeric(abs(op1$x-9.5)) < 10^-4, TRUE)
  
  op1 <- optim_DFP_NAPA(8, fun, xdat=c(8,11), tolerance=10^-9)
  expect_equal(as.numeric(abs(op1$x-9.5)) < 10^-6, TRUE)
  
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
    sum((par-centers)^4)
  }

  par     <- c(1,2,3,4,5)
  centers <- c(5,4,3,2,1)

  op1 <- optim_DFP_NAPA(starts = par, func = funND, centers=centers, tolerance = 10^-8, maxSteps = 200)
  expect_equal(as.numeric(abs(op1$x-centers)) < 10^-2, c(TRUE,TRUE,TRUE,TRUE,TRUE))
  
  if(require(redNMix)) {
    op1 <- redNMix::fit_red_Nmix_closed(nit = matrix(c(5,5,5),nrow=1), lambda_site_covariates = NULL, pdet_site_covariates = NULL, red = c(1), K = matrix(10, ncol=3), starts = c(1,0), method="DFP", tolerance=10^-4)
    op2 <- unmarked::pcount(formula = ~1 ~1, data = unmarked::unmarkedFramePCount(matrix(c(5,5,5),nrow=1)), K = 10, starts = c(1,0), se = FALSE)
    op3 <- redNMix::fit_red_Nmix_closed(nit = matrix(c(5,5,5),nrow=1), lambda_site_covariates = NULL, pdet_site_covariates = NULL, red = c(1), K = matrix(10, ncol=3), starts = c(1,0), method="BFGS", tolerance=10^-6)
    expect_equal(as.numeric(abs(exp(op1$x[1])-exp(op2@opt$par[1]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(exp(op1$x[1])-exp(op3$par[1]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(plogis(op1$x[2])-plogis(op2@opt$par[2]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(plogis(op1$x[2])-plogis(op3$par[2]))) < 10^-3, TRUE)
    
    op1 <- redNMix::fit_red_Nmix_closed(nit = matrix(c(5,6,4),nrow=1), lambda_site_covariates = NULL, pdet_site_covariates = NULL, red = c(1), K = matrix(10, ncol=3), starts = c(1,0), method="DFP", tolerance=10^-6, maxSteps=200)
    op2 <- unmarked::pcount(formula = ~1 ~1, data = unmarked::unmarkedFramePCount(matrix(c(5,6,4),nrow=1)), K = 10, starts = c(1,0), se = FALSE)
    op3 <- redNMix::fit_red_Nmix_closed(nit = matrix(c(5,6,4),nrow=1), lambda_site_covariates = NULL, pdet_site_covariates = NULL, red = c(1), K = matrix(10, ncol=3), starts = c(1,0), method="BFGS", tolerance=10^-8)
    expect_equal(as.numeric(abs(exp(op1$x[1])-exp(op2@opt$par[1]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(exp(op1$x[1])-exp(op3$par[1]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(plogis(op1$x[2])-plogis(op2@opt$par[2]))) < 10^-3, TRUE)
    expect_equal(as.numeric(abs(plogis(op1$x[2])-plogis(op3$par[2]))) < 10^-3, TRUE)
  }
  
})



test_that("optim_DFP_APA", {
  op1 <- optim_DFP_APA(15, func=function(x, precBits=64){(Rmpfr::mpfr(x,precBits)-10)^2})
  op2 <- optim(par = 15, fn = function(x){(x-10)^2}, method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-10, TRUE)
  
  op1 <- optim_DFP_APA(-15, func=function(x, precBits=64){(Rmpfr::mpfr(x,precBits)-10)^2})
  op2 <- optim(par = -15, fn = function(x){(x-10)^2}, method="BFGS")
  expect_equal(as.numeric(abs(op1$x-op2$par)) < 10^-10, TRUE)

  fun <- function(par, xdat, precBits=53) {
    l <- -1
    for(x in xdat) {
      l <- l * dpois_APA(x = x, lambda = par, precBits) 
    }
    return(l)
  }

  op1 <- optim_DFP_APA(starts = 10.0, func = fun, xdat=c(8,11), precBits=128, lineSearchMaxSteps = 50, keepValues = FALSE)
  expect_equal(as.numeric(abs(op1$x-mean(c(8,11)))) < 10^-6, TRUE)

  
  funND <- function(par, centers, precBits=53) {
    par     <- Rmpfr::mpfr(par, precBits)
    centers <- Rmpfr::mpfr(centers, precBits)

    sum((par-centers)^4)
  }

  precBits <- 53
  par     <- c(1,2,3,4,5)
  centers <- c(5,4,3,2,1)

  op1 <- optim_DFP_APA(starts = par, func = funND, centers=centers, precBits=precBits, maxSteps = 300, lineSearchMaxSteps = 200)
  expect_equal(as.numeric(abs(op1$x-centers)) < 10^-1, c(TRUE,TRUE,TRUE,TRUE,TRUE))
 
})
