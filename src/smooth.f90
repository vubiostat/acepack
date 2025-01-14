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

! This is an internal part of the smoother used by AVAS.
! See smothr.f90 and supsmu.90
!     
SUBROUTINE smooth (n,x,y,w,span,iper,vsmlsq,smo,acvr)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: x(n), y(n), w(n), span
  INTEGER,          INTENT(IN)    :: iper
  DOUBLE PRECISION, INTENT(IN)    :: vsmlsq
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), acvr(n)

  DOUBLE PRECISION xti,wt,fbw,xm,ym,var,cvar,fbo,tmp,xto,a,h,sy
  INTEGER          in,out,ibw, i,j,j0,it,jper
  
  xm   = 0.0
  ym   = 0.0
  var  = 0.0
  cvar = 0.0
  fbw  = 0.0
  jper = iabs(iper)
  ibw  = int(0.5*span*n+0.5)
  IF (ibw < 2) ibw=2
  it   = 2*ibw+1
  
  DO i=1,it
    j=i
    IF (jper==2) j=i-ibw-1
    xti=x(j)
    IF (j < 1) THEN
      j=n+j
      xti=x(j)-1.0
    END IF
    wt=w(j)
    fbo=fbw
    fbw=fbw+wt
    xm=(fbo*xm+wt*xti)/fbw
    ym=(fbo*ym+wt*y(j))/fbw
    tmp=0.0
    IF (fbo.gt.0.0) tmp=fbw*wt*(xti-xm)/fbo
    var=var+tmp*(xti-xm)
    cvar=cvar+tmp*(y(j)-ym)
  END DO
  
  DO j=1,n
    out = j-ibw-1
    in  = j+ibw
    IF (jper == 2 .or. (out >= 1 .and. in <= n)) THEN
      IF (out < 1) THEN
        out=n+out
        xto=x(out)-1.0
        xti=x(in)
      ELSE IF (in > n) THEN
        in=in-n
        xti=x(in)+1.0
        xto=x(out)
      ELSE
        xto=x(out)
        xti=x(in)
      END IF
 
      wt=w(out)
      fbo=fbw
      fbw=fbw-wt
      tmp=0.0
      IF (fbw > 0.0) tmp=fbo*wt*(xto-xm)/fbw
      var=var-tmp*(xto-xm)
      cvar=cvar-tmp*(y(out)-ym)
      xm=(fbo*xm-wt*xto)/fbw
      ym=(fbo*ym-wt*y(out))/fbw
      wt=w(in)
      fbo=fbw
      fbw=fbw+wt
      xm=(fbo*xm+wt*xti)/fbw
      ym=(fbo*ym+wt*y(in))/fbw
      tmp=0.0
      IF (fbo > 0.0) tmp=fbw*wt*(xti-xm)/fbo
      var=var+tmp*(xti-xm)
      cvar=cvar+tmp*(y(in)-ym)
    END IF
 
    a = 0.0
    IF (var > vsmlsq) a=cvar/var
    smo(j)=a*(x(j)-xm)+ym
    IF (iper > 0) THEN
      h=1.0/fbw
      IF (var > vsmlsq) h=h+(x(j)-xm)**2/var
      acvr(j)=abs(y(j)-smo(j))/(1.0-w(j)*h)
    END IF
  END DO
 
  j=1
  DO WHILE (j < n)
    j0=j
    sy=smo(j)*w(j)
    fbw=w(j)
    
    DO WHILE (j < n .and. x(j+1) <= x(j))
      j=j+1
      sy=sy+w(j)*smo(j)
      fbw=fbw+w(j)
    END DO
 
    IF (j > j0) THEN
      sy=sy/fbw
      smo(j0:j) = sy
    END IF

    j=j+1
  END DO

END SUBROUTINE smooth
