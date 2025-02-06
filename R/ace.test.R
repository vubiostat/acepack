  #############################################################################
 #
# This file is part of acepack.
#
# Copyright 2024,2025 Hajo Holzmann, Bernhard Klar
#
# Permission to use, copy, modify, distribute, and sell this software and
# its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and that
# both that copyright notice and this permission notice appear in
# supporting documentation. No representations are made about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
###############################################################################

#' @name ace.test
#' @title ACE permutation test of independence
#' @description Performs a permutation test of independence or association. The
#'   alternative hypothesis is that there is a x and y are dependent. 
#'   
#' Code authored by Hajo Holzmann, Bernhard Klar.
#' @param x a numeric vector, or a matrix or data frame with two columns.
#' @param y a vector with same length as x. Default is NULL.
#' @param nperm number of permutations. Default is 999.
#' @return a list containing the following:
#' \itemize{
#'   \item{\code{ace}}{ The value of the test statistic.}
#'   \item{\code{pval}}{ The *p*-value of the test.}
#' }
#' @export
#' @references
#' Holzmann, H., Klar, B. 2025. "Lancaster correlation - a new dependence measure
#' linked to maximum correlation". Scandinavian Journal of Statistics.
#' 52(1):145-169 <doi:10.1111/sjos.12733>
#' @importFrom arrangements npermutations
#' @importFrom arrangements permutations
#' 
#' @examples
#' 
#' n <- 200
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
#' cor.test(y[,1], y[,2], method = "spearman")
#' ace.test(y)
#' 
ace.test <- function(x, y = NULL, nperm = 999)
{ 
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  
  if (is.matrix(x))
  {
    y = x[,2]
    x = x[,1]
  }
  a       <- ace(x, y)
  ace.cor <- as.vector( cor(a$tx, a$ty) )
  n       <- length(x)

  ts <- ace.cor
  if (n <= 6)
  {
    nperm <- npermutations(n) #use all permutations
    perm  <- permutations(x)
    exact <- TRUE
  } else
  {
    perm  <- permutations(x, nsample = nperm)
    exact <- FALSE
  }
  tp <- rep(0, nperm)
  for (i in 1:nperm)
  {
    a     <- ace(perm[i,], y)
    tp[i] <- as.vector( cor(a$tx, a$ty) )
  }
  pval <- (sum(tp > ts) + 1) / (nperm + 1)
  
  structure(
    list(ace = ace.cor, pval = pval, exact=exact, n=nperm),
    class=c("ace.test", "list"))
}

#' @name summary.ace.test
#' @title ACE permutation test summary
#' @description A S3 function to produce a summary 
#'   of the results of an ace.test.
#' @param object the test to summarize
#' @param ... additional arguments (ignored)
#' @param digits Number of significant digits to round too.
#' @return a rounded ace.test object
#' @export
summary.ace.test <- function(object, ..., digits)
{
  object$ace  <- signif(object$ace,  digits)
  object$pval <- signif(object$pval, digits)
  
  object
}

#' @name summary.ace.test
#' @title ACE permutation test summary
#' @description A S3 function to produce a summary 
#'   of the results of an ace.test.
#' @param x the ace.test object to print
#' @param ... additional arguments to send to cat
#' @return NULL only outputs to CAT
#' @export
print.ace.test <- function(x, ...)
{
  if(x$exact)
  {
    cat("\nACE Exact Permutation Test of Independence\n", ...)
  } else {
    cat("\nACE Approximate Permutation Test of Independence\n", ...)
  }
  cat("\nalternative hypothesis: x and y are dependent\n", ...)
  cat("Ace correlation = ", x$ace, "\n", ...)
  pval <- format(x$pval, scientific=if(x$pval < 0.001) TRUE else FALSE)
  cat("p-value = ", pval, "\n", ...)
  cat("\n", ...)
}

