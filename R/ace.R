#' @name ace
#' @title Alternating Conditional Expectations
#'
#' @param x matrix; A matrix containing the independent variables.
#' @param y numeric; A vector containing the response variable.
#' @param wt numeric; An optional vector of weights.
#' @param cat integer; An optional integer vector specifying which variables
#'   assume categorical values.  Positive values in \code{cat} refer to columns
#'   of the \code{x} matrix and zero to the response variable.  Variables must
#'   be numeric, so a character variable should first be transformed with
#'   as.numeric() and then specified as categorical.
#' @param mon integer; An optional integer vector specifying which variables are
#'   to be transformed by monotone transformations.  Positive values
#'   in \code{mon} refer to columns of the \code{x} matrix and zero to the
#'   response variable.
#' @param lin integer; An optional integer vector specifying which variables are
#'   to be transformed by linear transformations. Positive values in \code{lin}
#'   refer to columns of the \code{x} matrix and zero to the response variable.
#' @param circ integer; An integer vector specifying which variables assume
#'   circular (periodic) values.  Positive values in \code{circ} refer to
#'   columns of the \code{x} matrix and zero to the response variable.
#' @param delrsq numeric(1); termination threshold. Iteration stops when
#'   R-squared changes by less than \code{delrsq} in 3 consecutive iterations
#'   (default 0.01).
#' @param control named list; control parameters to set. Documented at 
#' \code{\link{set_control}}.
#' @return
#'   A structure with the following components:
#'    \item{x}{the input x matrix.}
#'    \item{y}{the input y vector.}
#'    \item{tx}{the transformed x values.}
#'    \item{ty}{the transformed y values.}
#'    \item{rsq}{the multiple R-squared value for the transformed values.}
#'    \item{l}{the codes for cat, mon, ...}
#'
#' @description
#'    Uses the alternating conditional expectations algorithm to find the
#'    transformations of y and x that maximize the proportion of variation
#'    in y explained by x. When x is a matrix, it is transformed so that
#'   its columns are equally weighted when predicting y.
#'
#' @references
#'    Breiman and Friedman, Journal of the American Statistical
#'    Association (September, 1985).
#'
#'    The R code is adapted from S code for avas() by Tibshirani, in the
#'    Statlib S archive; the FORTRAN is a double-precision version of
#'   FORTRAN code by Friedman and Spector in the Statlib general archive.
#'
#' @examples
#' 
#' TWOPI <- 8*atan(1)
#' x <- runif(200,0,TWOPI)
#' y <- exp(sin(x)+rnorm(200)/2)
#' a <- ace(x,y)
#' par(mfrow=c(3,1))
#' plot(a$y,a$ty)  # view the response transformation
#' plot(a$x,a$tx)  # view the carrier transformation
#' plot(a$tx,a$ty) # examine the linearity of the fitted model
#' 
#' # example when x is a matrix
#' X1 <- 1:10
#' X2 <- X1^2
#' X <- cbind(X1,X2)
#' Y <- 3*X1+X2
#' a1 <- ace(X,Y)
#' plot(rowSums(a1$tx),a1$y)
#' (lm(a1$y ~ a1$tx)) # shows that the colums of X are equally weighted
#' 
#' # From D. Wang and M. Murphy (2005), Identifying nonlinear relationships
#' # regression using the ACE algorithm.  Journal of Applied Statistics,
#' # 32, 243-258.
#' X1 <- runif(100)*2-1
#' X2 <- runif(100)*2-1
#' X3 <- runif(100)*2-1
#' X4 <- runif(100)*2-1
#' 
#' # Original equation of Y:
#' Y <- log(4 + sin(3*X1) + abs(X2) + X3^2 + X4 + .1*rnorm(100))
#' 
#' # Transformed version so that Y, after transformation, is a
#' # linear function of transforms of the X variables:
#' # exp(Y) = 4 + sin(3*X1) + abs(X2) + X3^2 + X4
#' 
#' a1 <- ace(cbind(X1,X2,X3,X4),Y)
#' 
#' # For each variable, show its transform as a function of
#' # the original variable and the of the transform that created it,
#' # showing that the transform is recovered.
#' par(mfrow=c(2,1))
#' 
#' plot(X1,a1$tx[,1])
#' plot(sin(3*X1),a1$tx[,1])
#' 
#' plot(X2,a1$tx[,2])
#' plot(abs(X2),a1$tx[,2])
#' 
#' plot(X3,a1$tx[,3])
#' plot(X3^2,a1$tx[,3])
#' 
#' plot(X4,a1$tx[,4])
#' plot(X4,a1$tx[,4])
#' 
#' plot(Y,a1$ty)
#' plot(exp(Y),a1$ty)
#' 
#' @export
#' @useDynLib acepack, .registration=TRUE
ace <- function(x, y, wt = rep(1, nrow(x)), cat = NULL, mon = NULL, 
    lin = NULL, circ = NULL, delrsq = 0.01, control = NULL) 
{
  if(!is.null(control)) do.call(set_control, control)
  
  x  <- as.matrix(x)
  
  if (delrsq <= 0) stop("delrsq must be positive")

  iy <- ncol(x) + 1
  l  <- matrix(1, ncol = iy)
  
  if (!is.null(circ))
  {
    for (i in 1:length(circ))
    {
      if (circ[i] < 0 || circ[i] > ncol(x)) stop("bad circ= specification")

      nncol <- if (circ[i] == 0) iy else circ[i]

      if (l[nncol] != 2 & l[nncol] != 1) 
        stop("conflicting transformation specifications")

      l[nncol] <- 2
    }
  }
    
  if (length(mon) > 0)
  {
    for (i in 1:length(mon))
    {
      if (mon[i] < 0 || mon[i] > ncol(x)) stop("bad mon= specification")

      nncol <-  if (mon[i] == 0) iy else mon[i]

      if (l[nncol] != 3 && l[nncol] != 1) 
        stop("conflicting transformation specifications")

      l[nncol] <- 3
    }
  }
  
  if (length(lin)>0)
  {
    for (i in 1:length(lin))
    {
      if (lin[i] < 0 || lin[i] > ncol(x)) 
        stop("bad lin= specification")

      nncol <- if (lin[i] == 0) iy else lin[i]

      if (l[nncol] != 4 && l[nncol] != 1)
        stop("conflicting transformation specifications")
      
      l[nncol] <- 4
    }
  }
  
  if (length(cat))
  {
    for (i in 1:length(cat))
    {
      if (cat[i] < 0 || cat[i] > ncol(x)) stop("bad cat= specification")

      nncol <- if (cat[i] == 0) iy else cat[i]

      if (l[nncol] != 5 && l[nncol] != 1) 
        stop("conflicting transformation specifications")

      l[nncol] <- 5
    }
  }
  
  tx <- x
  ty <- y
  m  <- matrix(0, nrow = nrow(x), ncol = iy)
  z  <- matrix(0, nrow = nrow(x), ncol = 12)

  ns <- 1
  mode(x)      <- "double"
  mode(y)      <- "double"
  mode(tx)     <- "double"
  mode(ty)     <- "double"
  mode(wt)     <- "double"
  mode(delrsq) <- "double"
  mode(z)      <- "double"
  
  .Fortran("mace", p = as.integer(ncol(x)), n = as.integer(nrow(x)), 
    x = t(x), y = y, w = as.double(wt), l = as.integer(l), 
    delrsq = delrsq, ns = as.integer(ns), tx = tx, ty = ty, 
    rsq = double(1), ierr = integer(1), m = as.integer(m), 
    z = z, PACKAGE = "acepack")
}
