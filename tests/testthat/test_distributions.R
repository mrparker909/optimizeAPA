context("Testing Distribution Functions")

test_that("dpois_APA", {
  expect_equal(as.numeric(dpois_APA(x = 0:10, lambda = 5)),dpois(0:10, 5))
  expect_equal(as.numeric(dpois_APA(x = -1, lambda = 5)), 0)
  expect_equal(as.numeric(dpois_APA(x = -10:10, lambda = 5)),dpois(-10:10, 5))
  
  expect_equal(as.numeric(dpois_APA(x = seq(-1,1,0.5), lambda = 5)),dpois(seq(-1,1,0.5), 5))
})