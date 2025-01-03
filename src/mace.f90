  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 1985,2007 Jerome H. Friedman, Department of Statistics and
!   Stanford Linear Accelerator Center, Stanford University
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

! Estimate multiple optimal transformations for regression and
! correlation by alternating conditional expectation estimates.
!
! input:
!
!    n : number of observations.
!    p : number of predictor variables for each observation.
!    x(p,n) : predictor data matrix.
!    y(n) : response values for the observations.
!       missing values are signified by a value (response or
!       predictor) greater than or equal to big.
!       (see below - default big = 1.0e20)
!    w(n) : weights for the observations.
!    l(p+1) : flag for each variable.
!       l(1) through l(p) : predictor variables.
!       l(p+1) : response variable.
!       l(i)=0 => ith variable not to be used.
!       l(i)=1 => ith variable assumes orderable values.
!       l(i)=2 => ith variable assumes circular (periodic) values
!                 in the range (0.0,1.0) with period 1.0.
!       l(i)=3 => ith variable transformation is to be monotone.
!       l(i)=4 => ith variable transformation is to be linear.
!       l(i)=5 => ith variable assumes categorical (unorderable) values.
!   delrsq : termination threshold. iteration stops when
!       rsq changes less than delrsq in nterm
!       consecutive iterations (see below - default nterm=3).
!   ns : number of eigensolutions (sets of transformations).
!
! output:
!
!   tx(n,p,ns) : predictor transformations.
!      tx(j,i,k) = transformed value of ith predictor for jth obs
!                  for kth eigensolution.
!   ty(n,ns) = response transformations.
!      ty(j,k) = transformed response value for jth observation
!                for kth eigensolution.
!   rsq(ns) = fraction of variance(ty<y>)
!                       p
!         explained by sum tx(i)<x(i)>  for each eigensolution.
!                      i=1
!   ierr : error flag.
!      ierr = 0 : no errors detected.
!      ierr > 0 : error detected - see format statements below.
!
! scratch:
!
!    m(n,p+1), z(n,12) : internal working storage.
!
! Note: mace uses an iterative procedure for solving the optimization
!    problem. defa<starting transformations are ty(j,k)=y(j),
!    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns. other starting transformat
!    can be specified (if desired) for either the response and/or any of
!    the predictor variables. This is signaled by negating the
!    corresponding l(i) value and storing the starting transformed
!    values in the corresponding array (ty(j,k), tx(j,i,k)) before
!    calling mace.
!
SUBROUTINE mace (p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
  USE acedata
  IMPLICIT NONE

  INTEGER :: n,p,pp1,m(n,p+1),l(p+1)
  INTEGER :: ns,ierr,i,is,ism1,iter,j,js,k,nit,np,nt
  DOUBLE PRECISION :: rsqi
  DOUBLE PRECISION :: cmn, cmx
  DOUBLE PRECISION :: y(n),x(p,n),w(n),ty(n,ns),tx(n,p,ns)
  DOUBLE PRECISION :: z(n,12),ct(10),rsq(ns)
  DOUBLE PRECISION :: delrsq
  DOUBLE PRECISION :: sm,sv,sw,sw1
  
  ierr=0
  pp1=p+1
  sm=0.0
  sv=sm
  sw=sv
  sw1=sw
  
  ! Check if any element in the array l is out of the range [-5, 5]
  IF (any(l(1:pp1) < -5 .or. l(1:pp1) > 5)) THEN
    ierr = 6
    RETURN
  END IF
  
  IF (l(pp1) == 0) THEN
    ierr = 4
    RETURN
  END IF
  
  np = count(l /= 0)  ! Count the number of non-zero elements in l
  
  if (np <= 0) then
    ierr = 5
    return
  end if
  
  sw = sum(w)  ! Sum all the elements of array w
  
  if (sw .le. 0.0) then
    ierr = 1
    return
  end if
  
  ! Main Loop
  DO is=1,ns
  
    IF (l(pp1) > 0) THEN
      ty(1:n, is) = y(1:n)
    END IF
  
    DO i=1,p

      IF (l(i) == 0) EXIT
      
      IF (l(i) >= 0) THEN
        DO j=1,n
          tx(j,i,is)=x(i,j)
        END DO
      END IF
 
      DO j=1,n
        IF (tx(j,i,is) < big) THEN
          sm  = sm  + w(j)*tx(j,i,is)
          sw1 = sw1 + w(j)
        END IF
      END DO
      
      IF (sw1 > 0.0) EXIT
 
      DO j=1,n
        tx(j,i,is)=0.0
      END DO

      sm=0.0
      sw1=sm
      
    END DO
    
    IF (l(i) /= 0) THEN 
      sm=sm/sw1
   
      DO j=1,n
        IF (tx(j,i,is) >= big) THEN
          tx(j,i,is) = 0.0
        ELSE 
          tx(j,i,is) = tx(j,i,is)-sm
        END IF
      END DO
   
      sm=0.0
      sw1=sm
    END IF

    DO j=1,n
      tx(j,i,is)=0.0
    END DO
      
    DO j=1,n
      IF (ty(j,is) < big) THEN
        sm=sm+w(j)*ty(j,is)
        sw1=sw1+w(j)
      END IF
    END DO

    IF (sw1 <= 0.0) THEN
      ierr=1
      RETURN
    END IF
    
    sm = sm/sw1
    
    DO j=1,n
      IF (ty(j,is) >= big) THEN
        ty(j,is) = 0.0
      ELSE
        ty(j,is) = ty(j,is)-sm
      END IF
    END DO
    
    sv = (sv + sum(w*ty(:,is)**2))/sw
    
    IF (sv > 0.0) THEN
      sv = 1.0/dsqrt(sv)
    ELSE
      IF (l(pp1) <=0) THEN
        ierr = 3
      ELSE 
        ierr = 2
      ENDIF 
      RETURN
    ENDIF
    
    ty(:,is) = ty(:,is)*sv
 
    IF (is == 1) THEN
    
      DO j=1,n
        m(j,pp1)=j
        z(j,2)=y(j)
      END DO
 
      call sort (z(1,2),m(1,pp1),1,n)
      
      DO i=1,p
        IF (l(i) /=0) THEN
          DO j=1,n
            m(j,i)=j
            z(j,2)=x(i,j)
          END DO
          call sort (z(1,2),m(1,i),1,n)
        END IF
      END DO
  
    END IF
 
    call scail (p,n,w,sw,ty(1,is),tx(1,1,is),delrsq,p,z(1,5),z(1,6))
    
    rsq(is)=0.0
    iter=0
!      nterm=min0(nterm,10) ! FIXME Not allowed 
    nt=0
    
    DO i=1,min0(nterm,10)
      ct(i)=100.0
    END DO
    
    cmx = delrsq+1.0
    cmn = 0.0
    DO WHILE (cmx-cmn > delrsq .or. iter < maxit)
      iter=iter+1
      nit=0
      
 340  rsqi=rsq(is)
      nit=nit+1
      do 360 j=1,n
      z(j,5)=ty(j,is)
      do 350 i=1,p
      if (l(i).ne.0) z(j,5)=z(j,5)-tx(j,i,is)
 350  continue
 360  continue
      do 420 i=1,p
      if (l(i).eq.0) go to 420
      do 370 j=1,n
      k=m(j,i)
      z(j,1)=z(k,5)+tx(k,i,is)
      z(j,2)=x(i,k)
      z(j,4)=w(k)
 370  continue
      call smothr (iabs(l(i)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
      sm=0.0
      do 380 j=1,n
      sm=sm+z(j,4)*z(j,3)
 380  continue
      sm=sm/sw
      do 390 j=1,n
      z(j,3)=z(j,3)-sm
 390  continue
      sv=0.0
      do 400 j=1,n
      sv=sv+z(j,4)*(z(j,1)-z(j,3))**2
 400  continue
      sv=1.0-sv/sw
      if (sv.le.rsq(is)) go to 420
      rsq(is)=sv
      do 410 j=1,n
      k=m(j,i)
      tx(k,i,is)=z(j,3)
      z(k,5)=z(j,1)-z(j,3)
 410  continue
 420  continue
      if (np.eq.1.or.rsq(is)-rsqi.le.delrsq.or.nit.ge.maxit) go to 430
      go to 340
 430  do 450 j=1,n
      k=m(j,pp1)
      z(j,2)=y(k)
      z(j,4)=w(k)
      z(j,1)=0.0
      do 440 i=1,p
      if (l(i).ne.0) z(j,1)=z(j,1)+tx(k,i,is)
 440  continue
 450  continue
      call smothr (iabs(l(pp1)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
      if (is.le.1) go to 490
      ism1=is-1
      do 480 js=1,ism1
      sm=0.0
      do 460 j=1,n
      k=m(j,pp1)
      sm=sm+w(k)*z(j,3)*ty(k,js)
 460  continue
      sm=sm/sw
      do 470 j=1,n
      k=m(j,pp1)
      z(j,3)=z(j,3)-sm*ty(k,js)
 470  continue
 480  continue
 
 490  sm=0.0
      sv=sm
      DO j=1,n
        k=m(j,pp1)
        sm=sm+w(k)*z(j,3)
        z(k,2)=z(j,1)
      END DO
 
      sm=sm/sw
      
      DO j=1,n
        z(j,3)=z(j,3)-sm
        sv=sv+z(j,4)*z(j,3)**2
      END DO
      
      sv=sv/sw
      
      IF (sv <= 0.0) THEN
        ierr=3
        RETURN
      END IF

      sv=1.0/dsqrt(sv)
      
      DO j=1,n
        k=m(j,pp1)
        ty(k,is)=z(j,3)*sv
      END DO
 
      sv = sum(w(:)*(ty(:,is)-z(:,2))**2)
      rsq(is)=1.0-sv/sw
      nt=mod(nt,nterm)+1
      ct(nt)=rsq(is)
      DO i=1,nterm
        cmn=min(100.0, ct(i))
        cmx=max(-100.0,ct(i))
      END DO
    END DO

  END DO
  RETURN
END SUBROUTINE mace

