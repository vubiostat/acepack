context("ACE Transform")

set.seed(1) # For repeatability

x  <- matrix(runif(500)*2 - 1, ncol=5)
e  <- rnorm(100)
y  <- log(4 + sin(4*x[,1]) + abs(x[,2]) + x[,3]^2 + + x[,4]^3 + x[,5] + 0.1*e)
  


# D. Wang, Murphy M,. Estimating Optimal Transformations for Multiple
# Regression Using the ACE Algorithm.
# Journal of Data Science 2(2004), 329-346.
test_that("Estimates Multiple Transformations",
{
  y  <- log(4 + sin(4*x[,1]) + abs(x[,2]) + x[,3]^2 + x[,4]^3 + x[,5] + 0.1*e)

  expect_no_error(model <- ace(x, y))
  
  # Linear offset were computed using lm
  expect_true(max(sin(4*x[,1]) - model$tx[,1] - 0.003874) < 0.1)
  expect_true(max(abs(x[,2])   - model$tx[,2] - 0.481000) < 0.1)
  expect_true(max(x[,3]^2      - model$tx[,3] - 0.321443) < 0.1)
  expect_true(max(x[,4]^3      - model$tx[,4] - 0.039418) < 0.12)
  expect_true(max(x[,5]        - model$tx[,5] - 0.008231) < 0.1)
})

test_that("Estimates Multiple Transformations Specified via Formula",
{
  expect_no_error(model <- ace(y~x[,1]+x[,2]+x[,3]+x[,4]+x[,5]))
  
  # Linear offset were computed using lm
  expect_true(max(sin(4*x[,1]) - model$tx[,1] - 0.003874) < 0.1)
  expect_true(max(abs(x[,2])   - model$tx[,2] - 0.481000) < 0.1)
  expect_true(max(x[,3]^2      - model$tx[,3] - 0.321443) < 0.1)
  expect_true(max(x[,4]^3      - model$tx[,4] - 0.039418) < 0.12)
  expect_true(max(x[,5]        - model$tx[,5] - 0.008231) < 0.1)
  
  expect_no_error(
    model <- ace(y~x1+x2+x3+x4+x5, 
                 data.frame(x1=x[,1],
                            x2=x[,2],
                            x3=x[,3],
                            x4=x[,4],
                            x5=x[,5])))
  
  # Linear offset were computed using lm
  expect_true(max(sin(4*x[,1]) - model$tx[,1] - 0.003874) < 0.1)
  expect_true(max(abs(x[,2])   - model$tx[,2] - 0.481000) < 0.1)
  expect_true(max(x[,3]^2      - model$tx[,3] - 0.321443) < 0.1)
  expect_true(max(x[,4]^3      - model$tx[,4] - 0.039418) < 0.12)
  expect_true(max(x[,5]        - model$tx[,5] - 0.008231) < 0.1)
  
})


test_that("Handles categorical properly",
{
  y <- rnorm(100)
  x <- sample(1:3, 100, replace=TRUE)
  
  expect_no_error(result <- ace(x, y, cat=1))
  expect_equal(result$ierr, 0)
  
  expect_warning(ace(rep(1, 100), y, cat=1), "no variance")
  expect_warning(ace(x, rep(1, 100), cat=0), "no variance")
})

test_that("Will stop on error if specified",
{
  expect_error(ace(rep(1, 100), rnorm(100), cat=1, on.error=stop), "no variance") 
})