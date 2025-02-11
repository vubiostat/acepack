context("ACE Permutation Test (acetest)")

set.seed(1)

n <- 200
x <- matrix(rnorm(n*2), n)
nu <- 2
d200 <- x / sqrt(rchisq(n, nu)/nu) #multivariate t

set.seed(1)
n <- 5
x <- matrix(rnorm(n*2), n)
nu <- 2
d5 <- x / sqrt(rchisq(n, nu)/nu) #multivariate t

test_that("Repeats known results for n=200",
{
  expect_no_error(x <- acetest(d200))
  
  expect_equal(x$n, 999)
  expect_equal(x$ace, 0.5614263, tolerance=1e-6)
  expect_equal(x$pval, 0.001,    tolerance=1e-6)
})
  
test_that("Repeats known results for n=5",
{
  expect_no_error(x <- acetest(d5))
  
  expect_equal(x$n, 120)
  expect_equal(x$ace, 0.1809508,  tolerance=1e-6)
  expect_equal(x$pval, 0.8347107, tolerance=1e-6)
})

test_that("Pulls variable names",
{
  joe <- d5[,1]
  jim <- d5[,2]
  expect_no_error(x <- acetest(joe, jim))
  expect_equal(x$xname, "joe")
  expect_equal(x$yname, "jim")
})

test_that("Pulls matrix names",
{
  zed <- d5
  colnames(zed) <- c("sue", "bev")
  expect_no_error(x <- acetest(zed))
  expect_equal(x$xname, "sue")
  expect_equal(x$yname, "bev")
})

test_that("Pulls data.frame names",
{
  zed <- as.data.frame(d5)
  names(zed) <- c("sam", "pat")
  expect_no_error(x <- acetest(zed))
  expect_equal(x$xname, "sam")
  expect_equal(x$yname, "pat")
})

test_that("Accepts data.frame", 
{
  expect_no_error(x <- acetest(as.data.frame(d5)))
  
  expect_equal(x$n, 120)
  expect_equal(x$xname, 'V1')
})

test_that("Errors if matrix is not 2 columns.",
{
  expect_error(acetest(matrix(1:9, ncol=3)), 'must be 2 columns')
  expect_error(acetest(matrix(1:9, ncol=1)), 'must be 2 columns')
})

test_that("Cannot have a matrix 'x' and a y specified",
{
  expect_error(acetest(matrix(1:8, ncol=2), 1:10), 
    "Cannot have a matrix for 'x' and provide 'y'")
})
  
test_that("When x is not a matrix must have a y",
{
  expect_error(acetest(1:4), "Must supply both 'x' and 'y'")
})

test_that("'nperm' must be a positive integer",
{
  expect_error(acetest(1:4, 5:8, nperm=-1),  "'nperm' must be a positive integer")
  expect_error(acetest(1:4, 5:8, nperm='r'), "'nperm' must be a positive integer")
  expect_error(acetest(1:4, 5:8, nperm=2:3), "'nperm' must be a positive integer")
})

test_that("Length of 'x' and 'y' must be the same",
{
  expect_error(acetest(1:2, 5:10), "Length of 'x' and 'y' must be the same")
})

  