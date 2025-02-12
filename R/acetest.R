  #############################################################################
 #
# This file is part of acepack.
#
# Copyright 2024,2025 Hajo Holzmann, Bernhard Klar
# Copyright 2025 Shawn Garbett (edits and extensions)
#
# Permission to use, copy, modify, distribute, and sell this software and
# its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and that
# both that copyright notice and this permission notice appear in
# supporting documentation. No representations are made about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
###############################################################################

  
permn <- function(x, fun = NULL, ...)
{
# DATE WRITTEN: 23 Dec 1997          LAST REVISED:  23 Dec 1997
# AUTHOR:  Scott D. Chasalow (Scott.Chasalow@users.pv.wau.nl)
#
# DESCRIPTION:
#             Generates all permutations of the elements of x, in a minimal-
#	change order. If x is a	positive integer,  returns all permutations
#	of the elements of seq(x). If argument "fun" is not null,  applies
#	a function given by the argument to each point. "..." are passed
#	unchanged to the function given by argument fun, if any.
#
#	Returns a list; each component is either a permutation, or the
#	results of applying fun to a permutation.
#
# REFERENCE:
#	Reingold, E.M., Nievergelt, J., Deo, N. (1977) Combinatorial
#	Algorithms: Theory and Practice. NJ: Prentice-Hall. pg. 170.
#
# SEE ALSO:
#	sample, fact, combn, hcube, xsimplex
#
# EXAMPLE:
#	# Convert output to a matrix of dim c(6, 720)
#	t(array(unlist(permn(6)), dim = c(6, gamma(7))))
#
#	# A check that every element occurs the same number of times in each
#	# position
#	apply(t(array(unlist(permn(6)), dim = c(6, gamma(7)))), 2, tabulate, 
#		nbins = 6)
#
#	# Apply, on the fly, the diff function to every permutation
#	t(array(unlist(permn(6, diff)), dim = c(5, gamma(7))))
#
	if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) x <- seq(
			x)
	n <- length(x)
	nofun <- is.null(fun)
	out <- vector("list", gamma(n + 1))
	p <- ip <- seqn <- 1:n
	d <- rep(-1, n)
	d[1] <- 0
	m <- n + 1
	p <- c(m, p, m)
	i <- 1
	use <-  - c(1, n + 2)
	while(m != 1) {
		out[[i]] <- if(nofun) x[p[use]] else fun(x[p[use]], ...)
		i <- i + 1
		m <- n
		chk <- (p[ip + d + 1] > seqn)
		m <- max(seqn[!chk])
		if(m < n)
			d[(m + 1):n] <-  - d[(m + 1):n]
		index1 <- ip[m] + 1
		index2 <- p[index1] <- p[index1 + d[m]]
		p[index1 + d[m]] <- m
		tmp <- ip[index2]
		ip[index2] <- ip[m]
		ip[m] <- tmp
	}
	out
}

#' @name acetest
#' @title ACE permutation test of independence
#' @description Performs a permutation test of independence or association. The
#'   alternative hypothesis is that x and y are dependent. 
#'   
#' Code authored by Bernhard Klar, Shawn Garbett.
#' @param x a numeric vector, or a matrix or data frame with two columns. The
#'   first column is the 'y' and the second column is the 'x' when
#'   calling \code{\link{ace}}.
#' @param y a vector with same length as x. Default is NULL.
#' @param nperm number of permutations. Default is 999.
#' @param object S3 object of test results to dispatch.
#' @param digits Number of significant digits to round for summary.
#' @param acol for plot; color of the point estimate of correlation
#' @param xlim for plot;xlimit of histogram
#' @param col for plot;color of histogram bars
#' @param breaks for plot;number of breaks. Default to 100.
#' @param main for plot; main title of plot
#' @param xlab for plot; x-axis label
#' @param lwd for plot; line width of point estimate
#' @param ... additional arguments to pass to \code{cor}.
#' @seealso \code{\link{cor}}
#' @return a list containing the following:
#' \itemize{
#'   \item{\code{ace}} The value of the test statistic.
#'   \item{\code{pval}} The *p*-value of the test.
#' }
#' @references
#' Holzmann, H., Klar, B. 2025. "Lancaster correlation - a new dependence measure
#' linked to maximum correlation". Scandinavian Journal of Statistics.
#' 52(1):145-169 <doi:10.1111/sjos.12733>
#' @importFrom stats cor
#' @export
#' @rdname ace.test
#' @examples
#' 
#' n <- 200
#' z <- matrix(rnorm(2*n), n) / sqrt(rchisq(n, 2)/2)
#' x <- z[,2]; y <- z[,1]
#' cor.test(x, y, method="spearman")
#' acetest(x, y)
#' 
#' plot(acetest(z))
acetest <- function(x, y = NULL, nperm = 999, ...)
{ 
  if(is.data.frame(x)) x <- as.matrix(x)
  
  # Check user supplied parameters
  if (is.matrix(x) )
  {
    if (ncol(x) != 2) stop("Matrix 'x' must be 2 columns.")
    if (!is.null(y))    stop("Cannot have a matrix for 'x' and provide 'y'.")
  } else # x is not a matrix
  { 
    if (is.null(y))     stop("Must supply both 'x' and 'y' or a 2 column matrix 'x'.")
  } 
  if (!is.numeric(nperm) || nperm[1] <= 0 || length(nperm) != 1)
     stop("'nperm' must be a positive integer.")

  if (!is.null(y) && length(x) != length(y))
    stop("Length of 'x' and 'y' must be the same.")
  
  # Extract variable names
  xname <- as.character(substitute(x))
  yname <- as.character(substitute(y))
  if (is.matrix(x))
  {
    nm <- colnames(x)
    y = x[,2]
    x = x[,1]
    if(!is.null(nm))
    {
      xname <- nm[1]
      yname <- nm[2]
    }
  }

  if(is.null(yname) || identical(yname, character(0)) || yname == '') yname <- 'y'
  
  # Do the alternative hypothesis estimate
  a       <- ace(x, y)
  ace.cor <- as.vector( cor(a$tx, a$ty, ...) )
  n       <- factorial(length(x))

  if (n <= nperm) # Use all permutations
  {
    nperm <- n
    perm  <- permn(x)
    exact <- TRUE
    tp    <- vapply(
      1:nperm,
      function(i)
      {
        a <- ace(perm[[i]],y)
        cor(a$tx[,1], a$ty, ...)
      },
      numeric(1))
  } else # Only do a bootstrap approximation
  {
    exact <- FALSE
    tp    <- sapply(1:nperm, function(i) {
      a <- ace(sample(x),y)
      cor(a$tx[,1], a$ty, ...)
    })
  }
  pval <- (sum(tp > ace.cor) + 1) / (nperm + 1)
  
  structure(
    list(ace = ace.cor, pval = pval, exact=exact, n=nperm,
         tp  = tp, xname=xname, yname=yname),
    class=c("acetest", "list"))
}

#' @rdname ace.test
#' @export
summary.acetest <- function(object, ..., digits)
{
  object$ace  <- signif(object$ace,  digits)
  object$pval <- signif(object$pval, digits)
  
  object
}

#' @rdname ace.test
#' @export
print.acetest <- function(x, ...)
{
  if(x$exact)
  {
    cat("\nACE Exact Permutation Test of Independence\n", ...)
  } else {
    cat("\nACE Approximate Permutation Test of Independence\n", ...)
  }
  cat('\nalternative hypothesis:',x$xname,'and', x$yname, 'are dependent\n', ...)
  cat("Ace correlation \u03c1 =", x$ace, "\n", ...)
  pval <- format(x$pval, scientific=if(x$pval < 0.0001) TRUE else FALSE)
  if(1/(x$n+1) == x$pval)
  {
    cat("p-value <", pval, "\n", ...)
  } else
  {
    cat("p-value =", pval, "\n", ...)
  }
  cat("\n", ...)
  invisible(x)
}

#' @importFrom graphics hist
#' @importFrom graphics abline
#' @export
#' @rdname ace.test
plot.acetest <- function(
    x, 
    acol='blue',
    xlim=c(min(x$tp),max(c(x$tp, ceiling(x$ace*10)/10))),
    col='black',
    breaks=100,
    main='ACE Correlation Permutations',
    xlab=bquote(rho(.(x$xname),.(x$yname))),
    lwd=2,
    ...)
{
  hist(x$tp, xlim=xlim, col=col, breaks=breaks, main=main, xlab=xlab, ...)
  abline(v=x$ace, col=acol, lwd=lwd)
  invisible()
}