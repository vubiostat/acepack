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

SUBROUTINE ctsub(n,u,v,y,ty)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: n
  DOUBLE PRECISION, INTENT(IN)  :: u(n), v(n), y(n)
  DOUBLE PRECISION, INTENT(OUT) :: ty(n)

  INTEGER :: i,j
    
  DO i=1,n
    IF (y(i) <= u(1)) THEN
      ty(i)=(y(i)-u(1))*v(1)
    ELSE
      ty(i)=0.0
      DO j=1,n
        IF (y(i) <= u(j)) EXIT
        IF (j > 1) ty(i)=ty(i)+(u(j)-u(j-1))*(v(j)+v(j-1))/2
      END DO

      IF (y(i) <= u(n)) THEN
        ty(i)=ty(i)+.5*(y(i)-u(j-1))*(2*v(j-1)+(y(i)-u(j-1))*(v(j)-v(j-1))/(u(j)-u(j-1)))
      ELSE
        ty(i)=ty(i)+(y(i)-u(n))*v(n)
      END IF

    END IF
  END DO
END SUBROUTINE ctsub