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

SUBROUTINE bakfit(iter,delrsq,rsq,sw,l,z,m,x,ty,tx,w,n,p,np)
  USE acedata
  IMPLICIT NONE

  INTEGER,          INTENT(IN)    :: iter
  DOUBLE PRECISION, INTENT(IN)    :: delrsq
  DOUBLE PRECISION, INTENT(INOUT) :: rsq
  DOUBLE PRECISION, INTENT(IN)    :: sw
  INTEGER,          INTENT(IN)    :: l(p)
  DOUBLE PRECISION, INTENT(INOUT) :: z(n,17)
  INTEGER,          INTENT(IN)    :: m(n,p)
  DOUBLE PRECISION, INTENT(IN)    :: x(n,p)
  DOUBLE PRECISION, INTENT(OUT)   :: ty(n)
  DOUBLE PRECISION, INTENT(INOUT) :: tx(n,p)
  DOUBLE PRECISION, INTENT(IN)    :: w(n)
  INTEGER,          INTENT(IN)    :: n,p,np
  
  DOUBLE PRECISION :: sv, sm, rsqi
  INTEGER          :: nit, i, j, k

  CALL calcmu(n,p,l,z,tx)
  
  ty(:) = ty(:) - z(:,10)
  
  nit=0
  DO
    rsqi = rsq
    nit = nit+1
    
    DO i = 1,p
      IF (l(i).gt.0) THEN
        DO j = 1,n
          k = m(j,i)
          z(j,1) = ty(k)+tx(k,i)
          z(j,2) = x(k,i)
          z(j,7) = w(k)
        END DO
        
        CALL smothr(l(i),n,z(1,2),z,z(1,7),z(1,6),z(1,11))
        
        sm = 0.0
        DO j = 1,n
          sm = sm+z(j,7)*z(j,6)
        END DO
        sm = sm/sw
        DO j = 1,n
          z(j,6) = z(j,6)-sm
        END DO
        sv = 0.0
        DO j = 1,n
          sv = sv+z(j,7)*(z(j,1)-z(j,6))**2
        END DO
        sv = 1.0-sv/sw
        rsq = sv
        DO j = 1,n
          k = m(j,i)
          tx(k,i) = z(j,6)
          ty(k) = z(j,1)-z(j,6)
        END DO
      END IF
    END DO
    
    IF (np == 1 .or. abs(rsq-rsqi) <= delrsq .or. nit >= maxit) EXIT
  END DO
  
  IF (rsq == 0.0 .and. iter == 0) THEN
    DO i = 1,p
      IF (l(i) > 0) tx(:,i) = x(:,i)
    END DO
  END IF

END SUBROUTINE bakfit

