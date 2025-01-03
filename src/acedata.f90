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
! These procedure parameters can be changed in the calling routine
! by defining the above labeled common and resetting the values with
! executable statements.
!
! maxit : maximum number of iterations.
! nterm : number of consecutive iterations for which
!         rsq must change less than delcor for convergence.
! span, alpha : super smoother parameters (see below).
! big : a large representable floating point number.
!
MODULE acedata

  IMPLICIT NONE

  ! Define PARAMETERs for clarity using DOUBLE PRECISION
  DOUBLE PRECISION :: alpha = 0.0d0
  DOUBLE PRECISION :: big   = 1.0d20
  DOUBLE PRECISION :: span  = 0.0d0

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
