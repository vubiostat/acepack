\name{ace}
\alias{ace}
\title{Alternating Conditional Expectations}
\usage{
ace(x, y, wt = rep(1, nrow(x)), cat = NULL, mon = NULL, lin = NULL,
   circ = NULL, delrsq = 0.01)
}
\arguments{
    \item{x}{a matrix containing the independent variables.}
    \item{y}{a vector containing the response variable.}
    \item{wt}{an optional vector of weights.}
    \item{cat}{an optional integer vector specifying which variables
	assume categorical values.  Positive values in \code{cat} refer
	to columns of the \code{x} matrix and zero to the response
	variable.  Variables must be numeric, so a character variable
	should first be transformed with as.numeric() and then specified
	as categorical.}
    \item{mon}{an optional integer vector specifying which variables are
	to be transformed by monotone transformations.  Positive values
	in \code{mon} refer to columns of the \code{x} matrix and zero
	to the response variable.}
    \item{lin}{an optional integer vector specifying which variables are
	to be transformed by linear transformations.  Positive values in
	\code{lin} refer to columns of the \code{x} matrix and zero to
	the response variable.}
    \item{circ}{an integer vector specifying which variables assume
	circular (periodic) values.  Positive values in \code{circ}
	refer to columns of the \code{x} matrix and zero to the response
	variable.}
    \item{delrsq}{termination threshold.  Iteration stops when R-squared
	changes by less than \code{delrsq} in 3 consecutive iterations
	(default 0.01).}
}
\description{
    Uses the alternating conditional expectations algorithm to find the
    transformations of y and x that maximise the proportion of variation
    in y explained by x. When x is a matrix, it is transformed so that
    its columns are equally weighted when predicting y.
}
\value{
    A structure with the following components:
    \item{x}{the input x matrix.}
    \item{y}{the input y vector.}
    \item{tx}{the transformed x values.}
    \item{ty}{the transformed y values.}
    \item{rsq}{the multiple R-squared value for the transformed values.}
    \item{l}{the codes for cat, mon, ...}
    \item{m}{not used in this version of ace}
}
\references{
    Breiman and Friedman, Journal of the American Statistical
    Association (September, 1985).

    The R code is adapted from S code for avas() by Tibshirani, in the
    Statlib S archive; the FORTRAN is a double-precision version of
    FORTRAN code by Friedman and Spector in the Statlib general
    archive.
}
\examples{
TWOPI <- 8*atan(1)
x <- runif(200,0,TWOPI)
y <- exp(sin(x)+rnorm(200)/2)
a <- ace(x,y)
par(mfrow=c(3,1))
plot(a$y,a$ty)  # view the response transformation
plot(a$x,a$tx)  # view the carrier transformation
plot(a$tx,a$ty) # examine the linearity of the fitted model

# example when x is a matrix
X1 <- 1:10
X2 <- X1^2
X <- cbind(X1,X2)
Y <- 3*X1+X2
a1 <- ace(X,Y)
plot(rowSums(a1$tx),a1$y)
(lm(a1$y ~ a1$tx)) # shows that the colums of X are equally weighted

# From D. Wang and M. Murphy (2005), Identifying nonlinear relationships
# regression using the ACE algorithm.  Journal of Applied Statistics,
# 32, 243-258.
X1 <- runif(100)*2-1
X2 <- runif(100)*2-1
X3 <- runif(100)*2-1
X4 <- runif(100)*2-1

# Original equation of Y:
Y <- log(4 + sin(3*X1) + abs(X2) + X3^2 + X4 + .1*rnorm(100))

# Transformed version so that Y, after transformation, is a
# linear function of transforms of the X variables:
# exp(Y) = 4 + sin(3*X1) + abs(X2) + X3^2 + X4

a1 <- ace(cbind(X1,X2,X3,X4),Y)

# For each variable, show its transform as a function of
# the original variable and the of the transform that created it,
# showing that the transform is recovered.
par(mfrow=c(2,1))

plot(X1,a1$tx[,1])
plot(sin(3*X1),a1$tx[,1])

plot(X2,a1$tx[,2])
plot(abs(X2),a1$tx[,2])

plot(X3,a1$tx[,3])
plot(X3^2,a1$tx[,3])

plot(X4,a1$tx[,4])
plot(X4,a1$tx[,4])

plot(Y,a1$ty)
plot(exp(Y),a1$ty)
}
\keyword{models}
