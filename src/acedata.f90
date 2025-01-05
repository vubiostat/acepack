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
!   ChatGPT assisted refactor. Changed to a module got rid of implicit. 
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
! span  : double(1); ACE; Span to use in smoothing represents the fraction of
!         observations in smoothing window. Automatic span selection is
!         performed if set to 0.0. Default is 0.0 (automatic).
!
!         For small samples (n < 40) or if there are substantial serial
!         correlations between obserations close in x - value, then
!         a prespecified fixed span smoother (span > 0) should be
!         used. Reasonable span values are 0.3 to 0.5.
!
! alpha : double(1); AVAS; Controls high frequency (small span) penality
!         used with automatic span selection (base tone control). An
!         alpha < 0.0 or alpha > 10.0 results in no effect. Default is 5.0.
!
! big   : double(1); ACE and AVAS; a large representable floating point number.
!         Default is 1.0e20.
!
! sml   : double(1); AVAS; a small number. should be set so that (sml)**(10.0)
!         does not cause floating point underflow Default is 1e-4.
!
! spans : double(3); AVAS span values for the three running linear smoothers.
!           * spans(1) : tweeter span.  Default is 0.05.
!           * spans(2) : midrange span. Default is 0.2.
!           * spans(3) : woofer span.   Default is 0.5.
!         Warning: These span values should be changed only with care.
!
! eps   : double(1); AVAS used to numerically stabilize slope calculations
!         for running linear fits.
!
MODULE acedata

  IMPLICIT NONE

  DOUBLE PRECISION :: alpha = 5.0
  DOUBLE PRECISION :: big   = 1.0e20
  DOUBLE PRECISION :: span  = 0.0
  DOUBLE PRECISION :: sml   = 1e-4
  DOUBLE PRECISION :: eps   = 1e-3
  DOUBLE PRECISION :: spans(3) = (/0.05, 0.2, 0.5/)

  INTEGER :: maxit = 20
  INTEGER :: nterm = 3

CONTAINS

  SUBROUTINE set_alpha(a)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: a
    alpha = a
  END SUBROUTINE set_alpha
  
  SUBROUTINE set_big(b)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: b
    big = b
  END SUBROUTINE set_big
  
  SUBROUTINE set_sml(s)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: s
    sml = s
  END SUBROUTINE set_sml
  
  SUBROUTINE set_eps(e)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: e
    eps = e
  END SUBROUTINE set_eps
  
  SUBROUTINE set_spans(sps)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: sps(3)
    spans = sps
  END SUBROUTINE set_spans
  
  SUBROUTINE set_span(s)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: s
    span = s
  END SUBROUTINE set_span
  
  SUBROUTINE set_maxit(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m
    maxit = m
  END SUBROUTINE set_maxit
  
  SUBROUTINE set_nterm(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    nterm = n
  END SUBROUTINE set_nterm
  
END MODULE acedata
