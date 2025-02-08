context("ACE Permutation Test")


test_that("Test repeats known results for n=200",
{
  set.seed(1)
  
  n <- 200
  x <- matrix(rnorm(n*2), n)
  nu <- 2
  y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t

  expect_no_error(x <- ace.test(y))
  
  expect_equal(x$n, 999)
  expect_equal(x$ace, 0.5614263, tolerance=1e-6)
  expect_equal(x$pval, 0.001,    tolerance=1e-6)
})
  
test_that("Test repeats known results for n=5",
{
  set.seed(1)
  n <- 5
  x <- matrix(rnorm(n*2), n)
  nu <- 2
  y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t

  expect_no_error(x <- ace.test(y))
  
  expect_equal(x$n, 120)
  expect_equal(x$ace, 0.1809508,  tolerance=1e-6)
  expect_equal(x$pval, 0.8347107, tolerance=1e-6)
})