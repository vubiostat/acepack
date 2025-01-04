context("AVAS")

test_that("AVAS Creates finite output",
{
  set.seed(1)
  x <- runif(200,0,2*pi)
  y <- exp(sin(x)+rnorm(200)/2)
  a <- avas(x,y)
  
  expect_true(all(is.finite(a$tx)))
})

test_that("Estimates Multiple Transformations",
{
  set.seed(2) # For repeatability
  
  x  <- matrix(runif(500)*2 - 1, ncol=5)
  e  <- rnorm(100)
  
  y  <- log(4 + sin(4*x[,1]) + abs(x[,2]) + x[,3]^2 + + x[,4]^3 + x[,5] + 0.1*e)
  
  model <- avas(x, y)
  
  expect_equal(max(sin(4*x[,1]) - model$tx[,1]), 0.09213036, tol=1e-7)
  expect_equal(max(abs(x[,2])   - model$tx[,2]), 0.6359979,  tol=1e-7)
  expect_equal(max(x[,3]^2      - model$tx[,3]), 0.4453127,  tol=1e-7)
  expect_equal(max(x[,4]^3      - model$tx[,4]), 0.1625947,  tol=1e-7)
  expect_equal(max(x[,5]        - model$tx[,5]), 0.1094684,  tol=1e-7)
})