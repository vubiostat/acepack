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
! supporting documentation. no representations are made about the
! suitability of this software for any purpose.  it is provided "as is"
! without express or implied warranty.
!______________________________________________________________________________

SUBROUTINE calcmu(n,p,l,z,tx)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: n, p, l(p)
  DOUBLE PRECISION, INTENT(OUT) :: z(n, 17)
  DOUBLE PRECISION, INTENT(IN)  :: tx(n,p)
  
  INTEGER :: i

  z(:,10) = 0.0

  DO i=1,p
    IF (l(i) > 0) z(:,10) = z(:,10) + tx(:,i)
  END DO

END SUBROUTINE calcmu

