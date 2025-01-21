  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 1985,2007 Jerome H. Friedman
! Copyright 2016,2025 Shawn Garbett, Vanderbilt University Medical Center
!
! Permission to use, copy, modify, distribute, and sell this software and
! its documentation for any purpose is hereby granted without fee,
! provided that the above copyright notice appear in all copies and that
! both that copyright notice and this permission notice appear in
! supporting documentation. No representations are made about the
! suitability of this software for any purpose.  It is provided "as is"
! without express or implied warranty.
!______________________________________________________________________________

!
! 2025-01-02 Shawn Garbett
!   ChatGPT assisted refactor. Changed to a module. Elminated implicits.
!
! There are SET_* routines for each parameter available. 
! These parameters are used in the smoothing routines of ACE and AVAS. ACE and
! AVAS both have their own smoothing implementations. 
!
! maxit : integer(1); ACE and AVAS; Maximum number of iterations. Default is 20.
!
! nterm : integer(1); ACE and AVAS; Number of consecutive iterations for which
!         rsq must change less than delcor for convergence. Default is 3.
!
! span  : double(1); ACE and AVAS; Span to use in smoothing represents the
!         fraction of observations in smoothing window. Automatic span selection
!         is performed if set to 0.0. Default is 0.0 (automatic).
!
!         For small samples (n < 40) or if there are substantial serial
!         correlations between obserations close in x - value, then
!         a prespecified fixed span smoother (span > 0) should be
!         used. Reasonable span values are 0.3 to 0.5.
!
! spans : double(3); AVAS span values for the three running linear smoothers.
!           * spans(1) : tweeter span.  Default is 0.05.
!           * spans(2) : midrange span. Default is 0.2.
!           * spans(3) : woofer span.   Default is 0.5.
!         Warning: These span values should be changed only with great care.
!
! alpha : double(1); AVAS; Controls high frequency (small span) penality
!         used with automatic span selection (base tone control). An
!         alpha < 0.0 or alpha > 10.0 results in no effect. Default is 5.0.
!
! big   : double(1); ACE and AVAS; a large representable floating point number.
!         Default is 1.0e30.
!
! sml   : double(1); AVAS; A small number. Should be set so that (sml)**(10.0)
!         does not cause floating point underflow Default is 1e-30.
!
! eps   : double(1); AVAS; Used to numerically stabilize slope calculations
!         for running linear fits.
!
! References
!
! J Friedman, W Stuetzle. Smoothing of Scatterplots. Stanford Project Orion.
! July 1982
! 
! J Friedman. A Variable Span Smoother. LCS Technical Report No. 5. SLAC 
! PUB-3477. November 1984

MODULE acedata

  IMPLICIT NONE

  DOUBLE PRECISION :: alpha = 5.0
  DOUBLE PRECISION :: big   = 1.0e30
  DOUBLE PRECISION :: span  = 0.0
  DOUBLE PRECISION :: sml   = 1.0e-30  
  DOUBLE PRECISION :: eps   = 1e-3
  DOUBLE PRECISION :: spans(3) = (/0.05, 0.2, 0.5/)

  INTEGER :: maxit = 20
  INTEGER :: nterm = 3
  
  PUBLIC :: set_alpha
  PUBLIC :: set_big
  PUBLIC :: set_span
  PUBLIC :: set_sml
  PUBLIC :: set_eps
  PUBLIC :: set_spans

CONTAINS

  SUBROUTINE set_alpha(a)   BIND(C, name="set_alpha_") 
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: a
    alpha = a
  END SUBROUTINE set_alpha
  
  SUBROUTINE set_big(b)     BIND(C, name="set_big_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: b
    big = b
  END SUBROUTINE set_big
  
  SUBROUTINE set_sml(s)     BIND(C, name="set_sml_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: s
    sml = s
  END SUBROUTINE set_sml
  
  SUBROUTINE set_eps(e)     BIND(C, name="set_eps_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: e
    eps = e
  END SUBROUTINE set_eps
  
  SUBROUTINE set_spans(sps) BIND(C, name="set_spans_")
    USE, INTRINSIC :: iso_c_binding  
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: sps(3)
    spans = sps
  END SUBROUTINE set_spans
  
  SUBROUTINE set_span(s)    BIND(C, name="set_span_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: s
    span = s
  END SUBROUTINE set_span
  
  SUBROUTINE set_maxit(m)   BIND(C, name="set_maxit_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN) :: m
    maxit = m
  END SUBROUTINE set_maxit
  
  SUBROUTINE set_nterm(n)   BIND(C, name="set_nterm_")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN) :: n
    nterm = n
  END SUBROUTINE set_nterm
  
END MODULE acedata
