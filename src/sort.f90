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


! Puts into a the permutation vector which sorts v into
! increasing order. Only elements from ii to jj are considered.
! arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements.
!     
! This is a modification of CACM algorithm #347 which is a modified
! Hoare quicksort.
!
! Richard C. Singleton. (1969). Algorithm 347: an efficient algorithm for
! sorting with minimal storage [M1]. _Communications of the ACM_, 23(3), 185-7.
! doi: 10.1145/362875.362901.
  
SUBROUTINE sort (v,a,ii,jj)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT) :: v(*)
  INTEGER, INTENT(INOUT)          :: a(jj)
  INTEGER, INTENT(IN)             :: ii, jj

!      dimension a(jj),v(*)
      integer iu(20),il(20)
      integer t,tt,ij,j,k,l
      integer m,i
      double precision vt,vtt
      
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=a(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      a(k+1)=t
      v(k+1)=vt
      go to 100
      
END SUBROUTINE sort