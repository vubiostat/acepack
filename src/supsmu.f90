  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 1984,1985,2007 Jerome H. Friedman
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

! Inputs:
!     n : number of observations (x,y - pairs).
!     x(n) : ordered abscissa values.
!     y(n) : corresponding ordinate (response) values.
!     w(n) : weight for each (x,y) observation.
!     iper : periodic variable flag.
!       iper=1 => x is ordered interval variable.
!       iper=2 => x is a periodic variable with values
!                 in the range (0.0,1.0) and period 1.0.
! Outputs:
!     smo(n) : smoothed ordinate (response) values.
!     scratch:
!     scratch(n,7) : internal working storage.
! 
! See acedata for global control parameters
SUBROUTINE supsmu (n,x,y,w,iper,smo,scratch)
  USE acedata
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: x(n), y(n), w(n)
  INTEGER,          INTENT(IN)    :: iper
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), scratch(n,7)

  INTEGER          :: i, j, jper
  DOUBLE PRECISION :: h(1), sw, sy, a, scale, vsmlsq, resmin, f
  
  IF (x(n) <= x(1)) THEN
    sy = sum(w(:)*y(:))
    sw = sum(w(:))
   
    a=sy/sw
    smo(:) = a
    RETURN
  END IF
      
  i=n/4
  j=3*i
  scale=x(j)-x(i)
  DO WHILE (scale <= 0.0)
    IF (j < n) j=j+1
    IF (i > 1) i=i-1
    scale=x(j)-x(i)
  END DO
  
  vsmlsq=(eps*scale)**2
  jper=iper
  
  IF (iper==2 .and. (x(1) < 0.0 .or. x(n) > 1.0)) jper=1
  IF (jper < 1 .or. jper > 2) jper=1
  IF (span > 0.0) THEN
    CALL smooth (n,x,y,w,span,jper,vsmlsq,smo,scratch)
    RETURN
  END IF
  DO i=1,3
    CALL smooth (n,x,y,w,spans(i),jper,vsmlsq,scratch(1,2*i-1),scratch(1,7))
    CALL smooth (n,x,scratch(1,7),w,spans(2),-jper,vsmlsq,scratch(1,2*i),h)
  END DO
  DO j=1,n
    resmin=big
    DO i=1,3
      IF (scratch(j, 2*i) < resmin) THEN
        resmin=scratch(j,2*i)
        scratch(j,7)=spans(i)
      END IF
    END DO
    IF (alpha > 0.0 .and. alpha <= 10.0 .and. resmin < scratch(j,6)) THEN
      scratch(j,7) = scratch(j,7) + &
                (spans(3)-scratch(j,7))*max(sml,resmin/scratch(j,6))**(10.0-alpha)
    END IF
  END DO
  CALL smooth (n,x,scratch(1,7),w,spans(2),-jper,vsmlsq,scratch(1,2),h)
  DO j=1,n
    IF (scratch(j,2) <= spans(1)) scratch(j,2)=spans(1)
    IF (scratch(j,2) >= spans(3)) scratch(j,2)=spans(3)
    f=scratch(j,2)-spans(2)
    IF (f >= 0.0) THEN
      f=f/(spans(3)-spans(2))
      scratch(j,4)=(1.0-f)*scratch(j,3)+f*scratch(j,5)
    ELSE
      f=-f/(spans(2)-spans(1))
      scratch(j,4)=(1.0-f)*scratch(j,3)+f*scratch(j,1)
    END IF
  END DO
  CALL smooth (n,x,scratch(1,4),w,spans(1),-jper,vsmlsq,smo,h)
END SUBROUTINE supsmu
 