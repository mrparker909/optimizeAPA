context("Testing Gradient Algorithms")

test_that("grad_FD_NAPA", {
  expect_equal(grad_FD_NAPA(func = function(x){x^2}, x=0), 0, tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x^2}, x=1), 2, tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x^2}, x=2), 4, tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x^2}, x=-1), -2, tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x^2}, x=-2), -4, tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x[1]^2 + x[2]^2}, x=c(1,3)), c(2,6), tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x[1]^2 + x[2]^2}, x=c(-1,-3)), c(-2,-6), tolerance=10^-10)
  expect_equal(grad_FD_NAPA(func = function(x){x[1]^2 + x[2]^2}, x=c(-1,3)), c(-2,6), tolerance=10^-10)
})

test_that("grad_FD_APA", {
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x^2}, x=0)), 0, tolerance=10^-16)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x^2}, x=1)), 2, tolerance=10^-16)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x^2}, x=2)), 4, tolerance=10^-16)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x^2}, x=-1)), -2, tolerance=10^-16)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x^2}, x=-2)), -4, tolerance=10^-16)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x[1]^2 + x[2]^2}, x=c(1,3))), c(2,6), tolerance=10^-10)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x[1]^2 + x[2]^2}, x=c(-1,-3))), c(-2,-6), tolerance=10^-10)
  expect_equal(as.numeric(grad_FD_APA(func = function(x,precBits){x[1]^2 + x[2]^2}, x=c(-1,3))), c(-2,6), tolerance=10^-10)
  
})
