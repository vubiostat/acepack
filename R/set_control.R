  #############################################################################
 #
# This file is part of acepack.
#
# Copyright 1985,2007 Jerome H. Friedman
# Copyright 2016,2025 Shawn Garbett, Vanderbilt University Medical Center
#
# Permission to use, copy, modify, distribute, and sell this software and
# its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and that
# both that copyright notice and this permission notice appear in
# supporting documentation. No representations are made about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
###############################################################################

#' @name set_control
#' @title Set internal parameters that control ACE and AVAS algorithms
#' 
#' @param alpha numeric(1); AVAS; Controls high frequency (small span) penalty
#'    used with automatic span selection (base tone control). An
#'    alpha < 0.0 or alpha > 10.0 results in no effect. Default is 5.0.
#' @param big numeric(1); ACE and AVAS; a large floating point number.
#'    Default is 1.0e30.
#' @param sml numeric(1); AVAS; A small number. Should be set so that `(sml)**(10.0)`
#'    does not cause floating point underflow. Default is 1e-30.
#' @param span numeric(1); ACE and AVAS; Span to use in smoothing represents the
#'         fraction of observations in smoothing window. Automatic span selection
#'         is performed if set to 0.0. Default is 0.0 (automatic).
#'
#'         For small samples (n < 40) or if there are substantial serial
#'         correlations between observations close in x - value, then
#'         a specified fixed span smoother (span > 0) should be
#'         used. Reasonable span values are 0.3 to 0.5.
#' @param spans numeric(3); AVAS; span values for the three running linear smoothers.
#'  \describe{
#'    \item{"spans(1)"}{Tweeter span. Default is 0.05.}
#'    \item{"spans(2)"}{Midrange span. Default is 0.2.}
#'    \item{"spans(3)"}{Woofer span. Default is 0.5.}
#' }
#' Warning: These span values should be changed only with great care.
#' @param eps numeric(1); AVAS; Used to numerically stabilize slope calculations
#'         for running linear fits.
#' @param maxit integer(1); ACE and AVAS; Maximum number of iterations.
#'   Default is 20.
#' @param nterm integer(1); ACE and AVAS; Number of consecutive iterations for
#'   which rsq must change less than delcor for convergence. Default is 3.
#' @return NULL
#' @description
#' 
#' These parameters are used in the smoothing routines of ACE and AVAS. ACE and
#' AVAS both have their own smoothing implementations. This sets them globally
#' for the package. 
#' 
#' The default values are good for the vast majority of cases. This routine
#' is included to provide complete control to the user, but is rarely needed.
#' 
#' @examples
#' set_control(maxit=40)
#' set_control(maxit=20)
#' set_control(alpha=5.0)
#' set_control(big=1e30, sml=1e-30)
#' set_control(eps=1e-3)
#' set_control(span=0.0, spans=c(0.05, 0.2, 0.5))
#' set_control(maxit=20, nterm=3)
#' @export
set_control <- function(alpha=NULL, big  =NULL, span =NULL, sml  =NULL, 
                        eps  =NULL, spans=NULL, maxit=NULL, nterm=NULL)
{
  if(!is.null(alpha))
  {
    if(!inherits(alpha, "numeric") || length(alpha) != 1)
      stop("Misspecified alpha. Must be numeric(1).")
    if(alpha < 0.0 || alpha > 10.0)
      warning("Alpha outside of {0,10} skips use of parameter.")
    mode(alpha) <- "double"
    .Fortran("set_alpha", a=alpha, PACKAGE = "acepack")
  }
  if(!is.null(big))
  {
    if(!inherits(big, "numeric") || length(big) != 1)
      stop("Misspecified big. Must be numeric(1)")
    mode(big) <- "double"
    .Fortran("set_big",b=big, PACKAGE = "acepack")
  }
  if(!is.null(span))
  {
    if(!inherits(span, "numeric") || length(span) != 1)
      stop("Misspecified span. Must be numeric(1)")
    mode(span) <- "double"
    .Fortran("set_span",s=span, PACKAGE = "acepack")
  }
  if(!is.null(sml))
  {
    if(!inherits(sml, "numeric") || length(sml) != 1)
      stop("Misspecified sml. Must be numeric(1)")
    mode(sml) <- "double"
    .Fortran("set_sml",s=sml, PACKAGE = "acepack")
  }
  if(!is.null(eps))
  {
    if(!inherits(eps, "numeric") || length(eps) != 1)
      stop("Misspecified eps. Must be numeric(1)")
    mode(eps) <- "double"
    .Fortran("set_eps",e=eps, PACKAGE = "acepack")
  }
  if(!is.null(spans))
  {
    if(!inherits(spans, "numeric") || length(spans) != 3)
      stop("Misspecified spans. Must be numeric(3)")
    mode(spans) <- "double"
    .Fortran("set_spans",sps=spans, PACKAGE = "acepack")
  }
  if(!is.null(maxit))
  {
    if(!inherits(maxit, "numeric") || length(maxit) != 1)
      stop("Misspecified maxit. Must be numeric(1)")
    if(maxit <= 0)
      stop("Misspecified maxit. Must be larger than 0")
    .Fortran("set_maxit",m=as.integer(maxit), PACKAGE = "acepack")
  }
  if(!is.null(nterm))
  {
    if(!inherits(nterm, "numeric") || length(nterm) != 1)
      stop("Misspecified nterm. Must be numeric(1)")
    if(nterm <= 0)
      stop("Misspecified nterm. Must be larger than 0")
    .Fortran("set_nterm",m=as.integer(nterm), PACKAGE = "acepack")
  }
  
  invisible(NULL)
}