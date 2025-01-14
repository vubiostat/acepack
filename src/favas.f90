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

SUBROUTINE favas(p,n,x,y,w,l,delrsq,tx,ty,rsq,ierr,m,z,yspan,iter,iters)
  USE acedata
  IMPLICIT NONE
  
  INTEGER,          INTENT(IN)    :: p, n
  DOUBLE PRECISION, INTENT(IN)    :: x(n,p), y(n), w(n)
  INTEGER,          INTENT(IN)    :: l(p)
  DOUBLE PRECISION, INTENT(OUT)   :: delrsq, tx(n,p), ty(n), rsq
  INTEGER,          INTENT(OUT)   :: ierr, m(n, p)
  DOUBLE PRECISION, INTENT(OUT)   :: z(n, 17), yspan
  INTEGER,          INTENT(OUT)   :: iter
  DOUBLE PRECISION, INTENT(OUT)   :: iters(100,2)
  
  INTEGER          :: pp1, pp2, i, j, k, np, nt
  DOUBLE PRECISION :: sumlog, tres, rr, rnew, cmn, cmx
  DOUBLE PRECISION :: ct(10), rss, dof, sm, sv, sw, svx   
      
  ierr = 0
  pp1  = p + 1
  pp2  = p + 2
  np   = COUNT(l(1:p) > 0)
  
  sm = sum(w(:)*y(:))
  sv = sum(w(:)*y(:)**2)
  sw = sum(w(:))
  z(:,2) = y(:)
  m(1:n, pp1) = (/(j, j=1,n)/)
  sm = sm/sw
  sv = sv/sw-sm**2
  sv = 1.0/dsqrt(sv)
  z(:,1) = (y(:)-sm)*sv

  CALL sort(z(1,2), m(1,pp1), 1, n)
  
  DO i = 1,p
    IF (l(i) <= 0) CYCLE
    sm = sum(w(:)*x(:,i))/sw
    m(1:n, i) = (/(j, j=1,n)/)
    z(:,2) = x(:,i)
    CALL sort(z(1,2),m(1,i),1,n)
  END DO

  rsq = 0.0
  iter = 0
  nt = 0
  ct(:) = 100.0
  ty(:) = z(:, 1)
  z(:,9) = ty(:)
  CALL bakfit(iter,delrsq,rsq,sw,l,z,m,x,z(1,9),tx,w,n,p,np)
  
  sumlog=0
  
  DO
    iter = iter +1
    IF (l(pp1) /= 4) THEN
      CALL calcmu(n,p,l,z,tx)
      
      DO j=1,n
        tres=(ty(j)-z(j,10))
        IF (abs(tres) < 1e-10) tres=1e-10
        z(j,2)=log(sqrt(tres**2))
        m(j,pp2)=j
      END DO
  
      CALL sort(z(1,10),m(1,pp2),1,n)
      
      DO j=1,n
        k=m(j,pp2)
        z(j,4)=z(k,2)
        z(j,5)=w(k)
      END DO
      CALL rlsmo(z(1,10),z(1,4),z(1,5),yspan,dof,n,z(1,6),rss,z(1,7))
      
      DO j=1,n
        k=m(j,pp2)
        z(j,7)=exp(-z(j,6))
        sumlog=sumlog+n*(w(j)/sw)*2*z(j,6)
        z(j,8)=ty(k)
      END DO
  
      CALL ctsub(n,z(1,10),z(1,7),z(1,8),z(1,9))
      
      sm = sum(w(:)*z(:,9))
  
      DO j=1,n
        k=m(j,pp2)
        ty(k)=z(j,9)-sm/sw
      END DO
  
      sv  = sum((w(:)/sw)*ty(:)*ty(:))
      svx = sum((w(:)/sw)*z(:,10)*z(:,10))
  
      ty(:) = ty(:)/dsqrt(sv)
      
      DO j=1,n
        DO i=1,p
          IF (l(i) > 0) tx(j,i)=tx(j,i)/dsqrt(svx)
        END DO
      END DO
    END IF
 
    z(:,9) = ty(:)
    CALL bakfit(iter,delrsq,rsq,sw,l,z,m,x,z(1,9),tx,w,n,p,np)
    sumlog=sumlog+n*dlog(sv)
    CALL calcmu(n,p,l,z,tx)
    rr = sum((w(:)/sw)*(ty(:)-z(:,10))**2)
    rsq=1-rr
    rnew=sumlog+rr
    iters(iter,1)=iter
    iters(iter,2)=rsq
    nt = mod(nt,min0(nterm,10))+1
    ct(nt) = rsq
    cmn = MINVAL(ct(1:min0(nterm,10)))
    cmx = MAXVAL(ct(1:min0(nterm,10)))
    IF (cmx-cmn <= delrsq .or. iter >= maxit .or. l(pp1)==4) EXIT
  END DO
END SUBROUTINE favas