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

! This is the smoother used by AVAS
!
! J Friedman. A Variable Span Smoother. LCS Technical Report No. 5. SLAC 
! PUB-3477. November 1984
!
SUBROUTINE smothr (l, n, x, y, w, smo, scratch)
  USE acedata
  IMPLICIT NONE
  
  INTEGER,          INTENT(IN)    :: l, n
  DOUBLE PRECISION, INTENT(IN)    :: x(n), y(n), w(n)
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), scratch(n, 7)

  DOUBLE PRECISION :: sm, sw, a, b, d
  INTEGER          :: i, j, j0
      
  IF (l >= 5) THEN 
    j=1
    DO
      j0=j
      sm=w(j)*y(j)
      sw=w(j)
      DO WHILE (j < n .and. x(j+1) <= x(j))
        j=j+1
        sm=sm+w(j)*y(j)
        sw=sw+w(j)
      END DO
      sm=sm/sw
      smo(j0:j) = sm
      j=j+1
      IF (j > n) RETURN
    END DO
  END IF   
      
   IF (l == 4) THEN
    sm = sum(w(:)*x(:)*y(:))
    sw = sum(w(:)*x(:)**2)
    b  = sum(w(:)*x(:))
    d  = sum(w(:))
    a  = sm/(sw-(b**2)/d)
    b  = b/d
    smo(:) = a*(x(:)-b)
    RETURN
  END IF
      
  CALL supsmu (n,x,y,w,l,smo,scratch)
  
  IF (l /= 3) RETURN
  
  scratch(:,     1) = smo 
  scratch(n:-1:1,2) = smo
  CALL montne (scratch(1,1),n)
  CALL montne (scratch(1,2),n)
  
  sm = sum((smo(:)-scratch(:,1))**2)
  sw = sum((smo(:)-scratch(n:-1:1,2))**2)
 
  IF (sm < sw) THEN
    smo(:) = scratch(:, 1)
  ELSE
    smo(:) = scratch(n:-1:1,2)
  END IF
  j=1
  DO
    j0=j
    DO WHILE (j < n .and. smo(j+1)==smo(j)) 
      j=j+1
    END DO
 
    IF (j > j0) THEN
      a = 0.0
      IF (j0 > 1) a=0.5*(smo(j0)-smo(j0-1))
      b = 0.0
      IF (j < n) b=0.5*(smo(j+1)-smo(j))
      d = (a+b)/(j-j0)
      IF (a==0.0 .or. b==0.0) d=2.0*d
      IF (a==0.0) a=b
      DO i=j0,j
        smo(i)=smo(i)-a+d*(i-j0)
      END DO
    END IF
    
    j=j+1
    IF (j > n) EXIT
  END DO
  
  j=1
  DO 
    j0=j
    sm=smo(j)
    DO WHILE (j < n .and. x(j+1) <= x(j))
      j=j+1
      sm=sm+smo(j)
    END DO
    sm=sm/(j-j0+1)
    smo(j0:j) = sm
    j=j+1
    if (j.gt.n) RETURN
  END DO


END SUBROUTINE smothr
