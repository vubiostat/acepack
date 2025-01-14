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
  INTEGER,          INTENT(IN)    :: n,p,np
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
  
  DOUBLE PRECISION :: sv, sm, rsqi
  INTEGER          :: nit, i

  CALL calcmu(n,p,l,z,tx)
  
  ty(:) = ty(:) - z(:,10)
  
  nit=0
  DO
    rsqi = rsq
    nit  = nit+1
    
    DO i = 1,p
      IF (l(i) > 0) THEN
        z(:,1) = ty(m(:,i))+tx(m(:,i),i)
        z(:,2) = x(m(:,i),i)
        z(:,7) = w(m(:,i))
        CALL smothr(l(i),n,z(1,2),z,z(1,7),z(1,6),z(1,11))
        sm = sum(z(:,7)*z(:,6))/sw
        z(:,6) = z(:,6) - sm
        sv = 1.0-sum(z(:,7)*(z(:,1)-z(:,6))**2)/sw
        rsq = sv
        tx(m(:,i), i) = z(:,6)
        ty(m(:,i))    = z(:,1) - z(:,6)
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

