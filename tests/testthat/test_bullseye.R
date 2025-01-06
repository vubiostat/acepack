context("pair test from bullseye")

ace_cor <- function(x,y,handle.na=TRUE)
{
  if(handle.na){
    pick <- complete.cases(x, y)
    x <- x[pick]
    y <- y[pick]
  }
  cat <- NULL
  if (is.factor(x)){
    x <- as.numeric(x)
    cat <- 1
  }
  if (is.factor(y)) {
    y <- as.numeric(y)
    cat <- c(cat,0)
  }
  acepack::ace(x,y, cat=cat)
}

test_that("ace_cor works",
{
  expect_no_error(results <- ace_cor(iris$Sepal.Length, iris$Species))
  expect_equal(results$rsq, 0.7027773, tol=1e-6)
})


