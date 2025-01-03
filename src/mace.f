  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 2007 Jerome H. Friedman
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


      subroutine mace (p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
      use acedata
      implicit none
c
c   subroutine mace(p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
c------------------------------------------------------------------
c
c estimate multiple optimal transformations for regression and
c correlation by alternating conditional expectation estimates.
c
c version 3/28/85.
c
c breiman and friedman, journal of the american statistical
c association (september, 1985)
c
c coded  and copywrite (c) 1985 by:
c
c                        jerome h. friedman
c                     department of statistics
c                               and
c                stanford linear accelerator center
c                        stanford university
c
c all rights reserved.
c
c
c input:
c
c    n : number of observations.
c    p : number of predictor variables for each observation.
c    x(p,n) : predictor data matrix.
c    y(n) : response values for the observations.
c       missing values are signified by a value (response or
c       predictor) greater than or equal to big.
c       (see below - default, big = 1.0e20)
c    w(n) : weights for the observations.
c    l(p+1) : flag for each variable.
c       l(1) through l(p) : predictor variables.
c       l(p+1) : response variable.
c       l(i)=0 => ith variable not to be used.
c       l(i)=1 => ith variable assumes orderable values.
c       l(i)=2 => ith variable assumes circular (periodic) values
c                 in the range (0.0,1.0) with period 1.0.
c       l(i)=3 => ith variable transformation is to be monotone.
c       l(i)=4 => ith variable transformation is to be linear.
c       l(i)=5 => ith variable assumes categorical (unorderable) values.
c   delrsq : termination threshold. iteration stops when
c       rsq changes less than delrsq in nterm
c       consecutive iterations (see below - default, nterm=3).
c   ns : number of eigensolutions (sets of transformations).
c
c output:
c
c   tx(n,p,ns) : predictor transformations.
c      tx(j,i,k) = transformed value of ith predictor for jth obs
c                  for kth eigensolution.
c   ty(n,ns) = response transformations.
c      ty(j,k) = transformed response value for jth observation
c                for kth eigensolution.
c   rsq(ns) = fraction of variance(ty<y>)
c                       p
c         explained by sum tx(i)<x(i)>  for each eigensolution.
c                      i=1
c   ierr : error flag.
c      ierr = 0 : no errors detected.
c      ierr > 0 : error detected - see format statements below.
c
c scratch:
c
c    m(n,p+1), z(n,12) : internal working storage.
c
c note: mace uses an iterative procedure for solving the optimization
c    problem. default starting transformations are ty(j,k)=y(j),
c    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns. other starting transformat
c    can be specified (if desired) for either the response and/or any of
c    the predictor variables. this is signaled by negating the
c    corresponding l(i) value and storing the starting transformed
c    values in the corresponding array (ty(j,k), tx(j,i,k)) before
c    calling mace.
c
c------------------------------------------------------------------
c
      integer n,p,pp1,m(n,p+1),l(p+1)
      integer ns,ierr,i,is,ism1,iter,j,js,k,nit,np,nt
      double precision rsqi
      double precision cmn, cmx
      double precision y(n),x(p,n),w(n),ty(n,ns),tx(n,p,ns)
      double precision z(n,12),ct(10),rsq(ns)
      double precision delrsq
      double precision sm,sv,sw,sw1
      ierr=0
      pp1=p+1
      sm=0.0
      sv=sm
      sw=sv
      sw1=sw
      do 10 i=1,pp1
      if (l(i).ge.-5.and.l(i).le.5) go to 10
      ierr=6
 10   continue
      if (ierr.ne.0) return
      if (l(pp1).ne.0) go to 20
      ierr=4
      return
 20   np=0
      do 30 i=1,p
      if (l(i).ne.0) np=np+1
 30   continue
      if (np.gt.0) go to 40
      ierr=5
      return
 40   do 50 j=1,n
      sw=sw+w(j)
 50   continue
      if (sw.gt.0.0) go to 60
      ierr=1
      return
 60   do 580 is=1,ns
      do 70 j=1,n
      if (l(pp1).gt.0) ty(j,is)=y(j)
 70   continue
      do 170 i=1,p
      if (l(i).ne.0) go to 90
      do 80 j=1,n
      tx(j,i,is)=0.0
 80   continue
      go to 170
 90   if (l(i).le.0) go to 110
      do 100 j=1,n
      tx(j,i,is)=x(i,j)
 100  continue
 110  do 120 j=1,n
      if (tx(j,i,is).ge.big) go to 120
      sm=sm+w(j)*tx(j,i,is)
      sw1=sw1+w(j)
 120  continue
      if (sw1.gt.0.0) go to 140
      do 130 j=1,n
      tx(j,i,is)=0.0
 130  continue
      sm=0.0
      sw1=sm
      go to 170
 140  sm=sm/sw1
      do 160 j=1,n
      if (tx(j,i,is).ge.big) go to 150
      tx(j,i,is)=tx(j,i,is)-sm
      go to 160
 150  tx(j,i,is)=0.0
 160  continue
      sm=0.0
      sw1=sm
 170  continue
      do 180 j=1,n
      if (ty(j,is).ge.big) go to 180
      sm=sm+w(j)*ty(j,is)
      sw1=sw1+w(j)
 180  continue
      if (sw1.gt.0.0) go to 190
      ierr=1
      return
 190  sm=sm/sw1
      do 210 j=1,n
      if (ty(j,is).ge.big) go to 200
      ty(j,is)=ty(j,is)-sm
      go to 210
 200  ty(j,is)=0.0
 210  continue
      do 220 j=1,n
      sv=sv+w(j)*ty(j,is)**2
 220  continue
      sv=sv/sw
      if (sv.le.0.0) go to 230
      sv=1.0/dsqrt(sv)
      go to 260
 230  if (l(pp1).le.0) go to 240
      ierr=2
      go to 250
 240  ierr=3
 250  return
 260  do 270 j=1,n
      ty(j,is)=ty(j,is)*sv
 270  continue
      if (is.ne.1) go to 310
      do 280 j=1,n
      m(j,pp1)=j
      z(j,2)=y(j)
 280  continue
      call sort (z(1,2),m(1,pp1),1,n)
      do 300 i=1,p
      if (l(i).eq.0) go to 300
      do 290 j=1,n
      m(j,i)=j
      z(j,2)=x(i,j)
 290  continue
      call sort (z(1,2),m(1,i),1,n)
 300  continue
 310  call scail (p,n,w,sw,ty(1,is),tx(1,1,is),delrsq,p,z(1,5),z(1,6))
      rsq(is)=0.0
      iter=0
!      nterm=min0(nterm,10)
      nt=0
      do 320 i=1,nterm
      ct(i)=100.0
 320  continue
 330  iter=iter+1
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
      do 500 j=1,n
      k=m(j,pp1)
      sm=sm+w(k)*z(j,3)
      z(k,2)=z(j,1)
 500  continue
      sm=sm/sw
      do 510 j=1,n
      z(j,3)=z(j,3)-sm
      sv=sv+z(j,4)*z(j,3)**2
 510  continue
      sv=sv/sw
      if (sv.le.0.0) go to 520
      sv=1.0/dsqrt(sv)
      go to 530
 520  ierr=3
      return
 530  do 540 j=1,n
      k=m(j,pp1)
      ty(k,is)=z(j,3)*sv
 540  continue
      sv=0.0
      do 550 j=1,n
      sv=sv+w(j)*(ty(j,is)-z(j,2))**2
 550  continue
      rsq(is)=1.0-sv/sw
      nt=mod(nt,nterm)+1
      ct(nt)=rsq(is)
      cmn=100.0
      cmx=-100.0
      do 560 i=1,nterm
      cmn=min(cmn,ct(i))
      cmx=max(cmx,ct(i))
 560  continue
      if (cmx-cmn.le.delrsq.or.iter.ge.maxit) go to 570
      go to 330
 570  continue
 580  continue
      return
      end

