context("Testing lineSearch Algorithms")

test_that("lineSearch_NAPA movement", {
  Fu <- function(x) {x^2}
  x_curr <- 1
  ls <- lineSearch_NAPA(x_curr = x_curr, dk = -0.2, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(1)), -1)
  
  ls1 <- lineSearch_NAPA(x_curr = 2, dk = -0.1, lineSearchMaxSteps = 100, func = Fu)
  ls2 <- lineSearch_NAPA(x_curr = 1.9, dk = -0.1, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(ls2$f_next-ls1$f_next), -1)
  
  ls1 <- lineSearch_NAPA(x_curr = 2, dk = -3, lineSearchMaxSteps = 100, func = Fu)
  ls2 <- lineSearch_NAPA(x_curr = -1, dk = 3, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(ls2$f_next-ls1$f_next), -1)
  
  x_curr <- 1
  ls <- lineSearch_NAPA(x_curr = x_curr, dk = -grad_FD_NAPA(func = Fu, x_val = x_curr), lineSearchMaxSteps = 100, func = Fu)
  expect_equal(ls$x_next, 0, tolerance=10^-10)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(1)), -1)
  
  x_curr <- -1
  ls <- lineSearch_NAPA(x_curr = x_curr, dk = -grad_FD_NAPA(func = Fu, x_val = x_curr), lineSearchMaxSteps = 100, func = Fu)
  expect_equal(ls$x_next, 0, tolerance=10^-10)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(1)), -1)
  
  Fu <- function(x, center) {(x-center)^2}
  
  x_curr <- 5
  center <- 10
  ls <- lineSearch_NAPA(x_curr = x_curr, dk = -grad_FD_NAPA(func = Fu, x_val = x_curr, center=center), lineSearchMaxSteps = 100, func = Fu, center=center)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(x_curr, center)), -1)
  
  x_curr <- 15
  ls <- lineSearch_NAPA(x_curr = x_curr, dk = -grad_FD_NAPA(func = Fu, x_val = x_curr, center=center), lineSearchMaxSteps = 100, func = Fu, center=center)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(x_curr, center)), -1)
  
  Fu <- function(lamb, x) {
    -1*prod(dpois(x = x, lambda = lamb))
  }
  x_curr <- c(9,10,15)
  l_curr <- 10
  ls <- lineSearch_NAPA(x_curr = l_curr, dk = -grad_FD_NAPA(func = Fu, x_val = l_curr, x=x_curr), func = Fu, x=x_curr)
  expect_equal(sign(l_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(l_curr, x_curr)), -1)
  
  l_curr <- 12
  ls <- lineSearch_NAPA(x_curr = l_curr, dk = -grad_FD_NAPA(func = Fu, x_val = l_curr, x=x_curr), func = Fu, x=x_curr)
  expect_equal(sign(l_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(l_curr, x_curr)), -1)
})

test_that("lineSearch_APA movement", {
  
  Fu <- function(x, center=0, precBits=53) {(Rmpfr::mpfr(x,precBits)-center)^2}
  
  precBits <- 53
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, precBits = precBits), lineSearchMaxSteps = 100, func = Fu, precBits = precBits)
  expect_equal(as.numeric(ls$x_next), 0, tolerance=10^-8)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(x_curr)), -1)
  
  x_curr <- -1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, precBits = precBits), lineSearchMaxSteps = 100, func = Fu, precBits = precBits)
  expect_equal(as.numeric(ls$x_next), 0, tolerance=10^-8)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(x_curr)), -1)
  
  precBits <- 128
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, precBits = precBits), lineSearchMaxSteps = 100, func = Fu, precBits = precBits)
  expect_equal(as.numeric(ls$x_next), 0, tolerance=10^-10)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(x_curr)), -1)
  
  Fu <- function(x, center=5, precBits=53) {(Rmpfr::mpfr(x,precBits)-center)^4}
  
  precBits <- 53
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, precBits = precBits), lineSearchMaxSteps = 100, func = Fu, precBits = precBits)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(x_curr)), -1)
  
  precBits <- 53
  x_curr <- -1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, precBits = precBits), lineSearchMaxSteps = 100, func = Fu, precBits = precBits)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(x_curr)), -1)
  
  
  x_curr <- 5
  center <- 10
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, center=center, precBits=53), lineSearchMaxSteps = 100, func = Fu, center=center, precBits=53)
  expect_equal(sign(x_curr-ls$x_next), -1)
  expect_equal(sign(ls$f_next-Fu(x_curr, center)), -1)
  
  x_curr <- 15
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(func = Fu, x_val = x_curr, center=center, precBits=53), lineSearchMaxSteps = 100, func = Fu, center=center, precBits=53)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(x_curr, center)), -1)
  
  Fu <- function(x, precBits) {Rmpfr::mpfr(x, precBits=precBits)^2}
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -0.2, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(x_curr-ls$x_next), 1)
  expect_equal(sign(ls$f_next-Fu(1, 64)), -1)
  
  ls1 <- lineSearch_APA(x_curr = 2, dk = -0.1, lineSearchMaxSteps = 100, func = Fu)
  ls2 <- lineSearch_APA(x_curr = 1.9, dk = -0.1, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(ls2$f_next-ls1$f_next), -1)
  
  ls1 <- lineSearch_APA(x_curr = 2, dk = -3, lineSearchMaxSteps = 100, func = Fu)
  ls2 <- lineSearch_APA(x_curr = 0.5, dk = -3, lineSearchMaxSteps = 100, func = Fu)
  expect_equal(sign(ls2$f_next-ls1$f_next), -1)
  
  ls <- lineSearch_APA(x_curr = c(2.06), dk = c(-0.2), func=function(x, precBits) { 
    l <- -1
    for(xi in c(1,2,3)) {
      l <- l * dpois_APA(xi, x, precBits)
    }
    return(l)
    })
  expect_equal(sign(2.06-ls$x_next), 1)
  
  ls1 <- lineSearch_APA(x_curr = c(2.000025), dk = c(-0.2), func=function(x, precBits) { 
    l <- -1
    for(xi in c(1,2,3)) {
      l <- l * dpois_APA(xi, x, precBits)
    }
    return(l)
  })
  ls2 <- lineSearch_APA(x_curr = c(1.9999999999999999418), dk = c(10^-10), func=function(x, precBits) { 
    l <- -1
    for(xi in c(1,2,3)) {
      l <- l * dpois_APA(xi, x, precBits)
    }
    return(l)
  })
  expect_equal(sign(ls1$f_next-ls2$f_next), 1)
  
  
  Fu <- function(lamb, precBits=53) { 
    l <- -1
    for(xi in c(1,2,3)) {
      l <- l * dpois_APA(xi, lamb, precBits)
    }
    return(l)
  }
  precBits <- 53
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 3
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 30
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 2
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(abs(Fu(x_curr,precBits)-Fu(ls$x_next, precBits))<10^-8, TRUE)
  
  precBits <- 128
  x_curr <- 1
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 3
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 30
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(sign(Fu(x_curr,precBits)-Fu(ls$x_next, precBits)), 1)
  
  x_curr <- 2
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, precBits = precBits), func = Fu, precBits = precBits)
  expect_equal(abs(Fu(x_curr,precBits)-Fu(ls$x_next, precBits))<10^-8, TRUE)
  
  
  # simple 2D quadratic function optimization
  Fu <- function(par, centers, precBits=53) {
    centers <- Rmpfr::mpfr(centers, precBits)
    par <- Rmpfr::mpfr(par, precBits)
    centerx <- centers[1]
    centery <- centers[2]
    (par[1]-centerx)^2+(par[2]-centery)^2
  }
  precBits <- 53
  centers <- c(1,3)
  x_curr  <- c(1,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,3)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(-3,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(-3,1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  precBits <- 128
  centers <- c(1,3)
  x_curr  <- c(1,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,3)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(-3,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(-3,1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
  centers <- c(1,3)
  x_curr  <- c(3,-1)
  ls <- lineSearch_APA(x_curr = x_curr, dk = -grad_FD_APA(Fu, x_curr, centers=centers, precBits = precBits), func=Fu, centers=centers, precBits = precBits)
  expect_equal(sign(Fu(x_curr,centers,precBits)-Fu(ls$f_next,centers,precBits)), 1)
  
})