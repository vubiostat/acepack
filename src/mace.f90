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
!    problem. default starting transformations are ty(j,k)=y(j),
!    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns. other starting transformations
!    can be specified (if desired) for either the response and/or any of
!    the predictor variables. This is signaled by negating the
!    corresponding l(i) value and storing the starting transformed
!    values in the corresponding array (ty(j,k), tx(j,i,k)) before
!    calling mace.
!
SUBROUTINE mace (p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
  USE acedata
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN)           :: p, n
  DOUBLE PRECISION, INTENT(IN)  :: x(p,n),y(n),w(n)
  INTEGER, INTENT(IN)           :: l(p+1)
  DOUBLE PRECISION, INTENT(IN)  :: delrsq
  INTEGER, INTENT(IN)           :: ns
  ! Outputs
  DOUBLE PRECISION, INTENT(OUT) :: tx(n,p,ns), ty(n,ns)
  DOUBLE PRECISION, INTENT(OUT) :: rsq(ns)
  INTEGER, INTENT(OUT)          :: ierr
  ! Scratch provided
  INTEGER, INTENT(OUT)          :: m(n,p+1)
  DOUBLE PRECISION, INTENT(OUT) :: z(n,12)
  
  ! Internal Variables
  INTEGER          :: pp1, i,is,ism1,iter,j,js,k,nit,np,nt
  DOUBLE PRECISION :: rsqi
  DOUBLE PRECISION :: cmn, cmx
  DOUBLE PRECISION :: ct(10)
  DOUBLE PRECISION :: sm,sv,sw,sw1
  
  ierr=0
  pp1=p+1
  sm=0.0
  sv=0.0
  sw=0.0
  sw1=0.0
  
  IF (any(l(1:pp1) < -5 .or. l(1:pp1) > 5)) THEN
    ierr = 6
    RETURN
  END IF
  
  IF (l(pp1) == 0) THEN
    ierr=4
    RETURN
  END IF

  np = count( l /= 0)

  IF (np <= 0) THEN
    ierr=5
    RETURN
  END IF
  
  sw = sum(w(:))
  IF (sw <= 0.0) THEN
    ierr=1
    RETURN
  END IF
      
  DO is=1,ns
    IF (l(pp1) > 0) ty(:,is)=y(:)

    DO i=1,p
      IF (l(i) == 0) THEN
        tx(:, i, is)=0.0
        CYCLE
      END IF
      
      IF (l(i) > 0) THEN
        tx(:,i,is)=x(i,:)
      END IF
 
      DO j=1,n
        IF (tx(j,i,is) < big) THEN
          sm=sm+w(j)*tx(j,i,is)
          sw1=sw1+w(j)
        END IF
      END DO
 
      IF (sw1 <= 0.0) THEN
        tx(:,i,is) = 0.0
        sm=0.0
        sw1=sm
        CYCLE
      END IF
      
      sm=sm/sw1
      
      DO j=1,n
        IF (tx(j,i,is) < big) THEN 
          tx(j,i,is)=tx(j,i,is)-sm
        ELSE
          tx(j,i,is)=0.0
        END IF
      END DO
 
      sm=0.0
      sw1=sm
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
    
    sm=sm/sw1
    
    DO j=1,n
      IF (ty(j,is) < big) THEN
        ty(j,is)=ty(j,is)-sm
      ELSE 
        ty(j,is)=0.0
      END IF
    END DO
    
    sv = (sv + sum(w(:)*ty(:,is)**2))/sw

    IF (sv <= 0.0) THEN
      IF (l(pp1) <= 0) THEN
        ierr=3
      ELSE
        ierr=2
      END IF
      RETURN
    END IF
    
    sv = 1.0/dsqrt(sv)
    ty(:,is) = ty(:,is)*sv

    IF (IS == 1) THEN
      m(1:n, pp1) = reshape([(i, i=1,n)], shape=[n]) ! 1:n
      z(:, 2) = y(:)
      CALL sort (z(1,2),m(1,pp1),1,n)
      DO i=1,p
        IF (l(i) /= 0) THEN
          m(1:n,i) = reshape([(i, i=1,n)], shape=[n]) ! 1:n
          z(:,2) = x(i,:)
          CALL sort (z(1,2),m(1,i),1,n)
        END IF
      END DO
    END IF
    
    CALL scail (p,n,w,sw,ty(1,is),tx(1,1,is),delrsq,p,z(1,5),z(1,6))
    rsq(is)=0.0
    iter=0
    nt=0
    
    ct(1:min0(nterm,10)) = 100.0
    
    DO ! Until maxit or convergence
      iter=iter+1
      nit=0
      
      DO 
        rsqi=rsq(is)
        nit=nit+1
        DO j=1,n
          z(j,5)=ty(j,is)
          DO i=1,p
            if (l(i).ne.0) z(j,5)=z(j,5)-tx(j,i,is)
          END DO
        END DO
        
        DO i=1,p
          IF (l(i) == 0) CYCLE
          DO j=1,n
            k=m(j,i)
            z(j,1)=z(k,5)+tx(k,i,is)
            z(j,2)=x(i,k)
            z(j,4)=w(k)
          END DO
          CALL smothr (iabs(l(i)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
          sm = sum(z(:,4)*z(:,3))/sw
          z(:,3) = z(:,3)-sm
          sv = sum(z(:,4)*(z(:,1)-z(:,3))**2)
          sv=1.0-sv/sw
          IF (sv > rsq(is)) THEN
            rsq(is)=sv
            DO j=1,n
              k=m(j,i)
              tx(k,i,is)=z(j,3)
              z(k,5)=z(j,1)-z(j,3)
            END DO
          END IF
        END DO
        IF (np==1 .or. rsq(is)-rsqi <= delrsq .or. nit >= maxit) EXIT
      END DO
 
      DO j=1,n
        k=m(j,pp1)
        z(j,2)=y(k)
        z(j,4)=w(k)
        z(j,1)=0.0
        DO i=1,p
          if (l(i) /= 0) z(j,1)=z(j,1)+tx(k,i,is)
        END DO
      END DO
 
      CALL smothr (iabs(l(pp1)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
      IF (is > 1) THEN
        ism1 = is-1
        DO js=1,ism1
          sm = sum(w(m(:,pp1))*z(:,3)*ty(m(:,pp1),js))/sw
          z(:,3) = z(:,3)-sm*ty(m(:,pp1),js)
        END DO
      END IF

      sm=sum(w(m(:,pp1))*z(:,3)) / sw
      z(m(:,pp1),2)=z(:,1)
      z(:,3)=z(:,3)-sm
      sv=sum(z(:,4)*z(:,3)**2) / sw
      IF (sv <= 0.0) THEN
        ierr=3
        RETURN
      END IF
      sv=1.0/dsqrt(sv)
      ty(m(:,pp1),is)=z(:,3)*sv
      sv=sum(w(:)*(ty(:,is)-z(:,2))**2)
      rsq(is)=1.0-sv/sw
      nt=mod(nt,min0(nterm,10))+1
      ct(nt)=rsq(is)
      cmn = minval(ct(1:min0(nterm, 10)))
      cmx = maxval(ct(1:min0(nterm, 10)))
      IF (cmx-cmn <= delrsq .or. iter >= maxit) EXIT
    END DO
  END DO
END SUBROUTINE mace
