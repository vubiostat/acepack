context("AVAS")

test_that(
  "AVAS Creates finite output",
{
  x <- runif(200,0,2*pi)
  y <- exp(sin(x)+rnorm(200)/2)
  a <- avas(x,y)
  
  expect_true(all(is.finite(a$tx)))
})