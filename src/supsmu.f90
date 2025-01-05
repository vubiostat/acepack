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
! Note:
!     For small samples (n < 40) or if there are substantial serial
!     correlations between obserations close in x - value, then
!     a prespecified fixed span smoother (span > 0) should be
!     used. Reasonable span values are 0.3 to 0.5.

SUBROUTINE supsmu (n,x,y,w,iper,smo,scratch)
  USE acedata
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: x(n), y(n), w(n)
  INTEGER,          INTENT(IN)    :: iper
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), scratch(n,7)

  INTEGER :: i, j, jper
  DOUBLE PRECISION :: h(1), sw, sy, a, scale, vsmlsq, resmin, f
  
      if (x(n).gt.x(1)) go to 30
      sy=0.0
      sw=sy
      do 10 j=1,n
         sy=sy+w(j)*y(j)
         sw=sw+w(j)
 10   continue
      a=sy/sw
      do 20 j=1,n
         smo(j)=a
 20   continue
      return
 30   i=n/4
      j=3*i
      scale=x(j)-x(i)
 40   if (scale.gt.0.0) go to 50
      if (j.lt.n) j=j+1
      if (i.gt.1) i=i-1
      scale=x(j)-x(i)
      go to 40
 50   vsmlsq=(eps*scale)**2
      jper=iper
      if (iper.eq.2.and.(x(1).lt.0.0.or.x(n).gt.1.0)) jper=1
      if (jper.lt.1.or.jper.gt.2) jper=1
      if (span.le.0.0) go to 60
      call smooth (n,x,y,w,span,jper,vsmlsq,smo,scratch)
      return
 60   do 70 i=1,3
         call smooth (n,x,y,w,spans(i),jper,vsmlsq,scratch(1,2*i-1),scratch(1,7))
         call smooth (n,x,scratch(1,7),w,spans(2),-jper,vsmlsq,scratch(1,2*i),h)
 70   continue
      do 90 j=1,n
         resmin=big
         do 80 i=1,3
            if (scratch(j,2*i).ge.resmin) go to 80
            resmin=scratch(j,2*i)
            scratch(j,7)=spans(i)
 80      continue
         IF (alpha.gt.0.0.and.alpha.le.10.0.and.resmin.lt.scratch(j,6)) THEN
           scratch(j,7) = scratch(j,7) + &
                     (spans(3)-scratch(j,7))*max(sml,resmin/scratch(j,6))**(10.0-alpha)
         END IF
 90   continue
      call smooth (n,x,scratch(1,7),w,spans(2),-jper,vsmlsq,scratch(1,2),h)
      do 110 j=1,n
         if (scratch(j,2).le.spans(1)) scratch(j,2)=spans(1)
         if (scratch(j,2).ge.spans(3)) scratch(j,2)=spans(3)
         f=scratch(j,2)-spans(2)
         if (f.ge.0.0) go to 100
         f=-f/(spans(2)-spans(1))
         scratch(j,4)=(1.0-f)*scratch(j,3)+f*scratch(j,1)
         go to 110
 100     f=f/(spans(3)-spans(2))
         scratch(j,4)=(1.0-f)*scratch(j,3)+f*scratch(j,5)
 110  continue
      call smooth (n,x,scratch(1,4),w,spans(1),-jper,vsmlsq,smo,h)
      return
END SUBROUTINE supsmu
 