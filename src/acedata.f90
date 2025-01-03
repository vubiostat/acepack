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
  DOUBLE PRECISION, PARAMETER :: alpha = 0.0d0
  DOUBLE PRECISION, PARAMETER :: big   = 1.0d20
  DOUBLE PRECISION, PARAMETER :: span  = 0.0d0

  INTEGER, PARAMETER :: maxit = 20
  INTEGER, PARAMETER :: nterm = 3

END MODULE acedata
