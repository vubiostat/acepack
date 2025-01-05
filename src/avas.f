
      block data avasdata
      integer itape, maxit, nterm
      double precision span, alpha
      double precision big, sml, eps
      double precision spans(3)
      common /parms/ span,alpha,itape,maxit,nterm
      common /spans/ spans /consts/ big,sml,eps
c------------------------------------------------------------------
c
c these procedure parameters can be changed in the calling routine
c by defining the above labeled common and resetting the values with
c executable statements.
c
c itape : fortran file number for printer output.
c         (itape.le.0 => no printer output.)
c maxit : maximum number of iterations.
c nterm : number of consecutive iterations for which
c         rsq must change less than delcor for convergence.
c span, alpha : super smoother parameters.
c   (see - friedman and stuetzle, reference above.)
c
c------------------------------------------------------------------
      data itape,maxit,nterm,span,alpha /-6,20,3,0.0,5.0/
c---------------------------------------------------------------
c
c this sets the compile time (default) values for various
c internal parameters :
c
c spans : span values for the three running linear smoothers.
c spans(1) : tweeter span.
c spans(2) : midrange span.
c spans(3) : woofer span.
c (these span values should be changed only with care.)
c big : a large representable floating point number.
c sml : a small number. should be set so that (sml)**(10.0) does
c       not cause floating point underflow.
c eps : used to numerically stabilize slope calculations for
c       running linear fits.
c
c these parameter values can be changed by declaring the
c relevant labeled common in the main program and resetting
c them with executable statements.
c
c-----------------------------------------------------------------
      data spans,big,sml,eps /0.05,0.2,0.5,1.0e20,1.0e-4,1.0e-3/
      end



      subroutine supsmu (n,x,y,w,iper,span,alpha,smo,sc)
c------------------------------------------------------------------
c     
c     super smoother (friedman and stuetzle, 1984).
c     
c     version 3/10/84
c     
c     coded by: j. h. friedman
c     department of statistics and
c     stanford linear accelerator center
c     stanford university
c     stanford ca. 94305
c     
c     input:
c     n : number of observations (x,y - pairs).
c     x(n) : ordered abscissa values.
c     y(n) : corresponding ordinate (response) values.
c     w(n) : weight for each (x,y) observation.
c     iper : periodic variable flag.
c     iper=1 => x is ordered interval variable.
c     iper=2 => x is a periodic variable with values
c     in the range (0.0,1.0) and period 1.0.
c     span : smoother span (fraction of observations in window).
c     span=0.0 => automatic (variable) span selection.
c     alpha : controles high frequency (small span) penality
c     used with automatic span selection (base tone control).
c     (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
c     output:
c     smo(n) : smoothed ordinate (response) values.
c     scratch:
c     sc(n,7) : internal working storage.
c     
c     note:
c     for small samples (n < 40) or if there are substantial serial
c     correlations between obserations close in x - value, then
c     a prespecified fixed span smoother (span > 0) should be
c     used. reasonable span values are 0.3 to 0.5.
c     
c------------------------------------------------------------------
      implicit none
      integer n
      double precision x(n),y(n),w(n),smo(n),sc(n,7)
      double precision big,sml,eps
      double precision span,alpha
      integer i,j,jper,iper
      double precision spans(3)
      common /spans/ spans /consts/ big,sml,eps
      double precision h(1),sw,sy,a,scale,vsmlsq,resmin,f
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
      call smooth (n,x,y,w,span,jper,vsmlsq,smo,sc)
      return
 60   do 70 i=1,3
         call smooth (n,x,y,w,spans(i),jper,vsmlsq,sc(1,2*i-1),sc(1,7))
         call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2*i),h)
 70   continue
      do 90 j=1,n
         resmin=big
         do 80 i=1,3
            if (sc(j,2*i).ge.resmin) go to 80
            resmin=sc(j,2*i)
            sc(j,7)=spans(i)
 80      continue
         if (alpha.gt.0.0.and.alpha.le.10.0.and.resmin.lt.sc(j,6))
     1     sc(j,7) = sc(j,7) +
     2        (spans(3)-sc(j,7))*max(sml,resmin/sc(j,6))**
     3        (10.0-alpha)
 90   continue
      call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2),h)
      do 110 j=1,n
         if (sc(j,2).le.spans(1)) sc(j,2)=spans(1)
         if (sc(j,2).ge.spans(3)) sc(j,2)=spans(3)
         f=sc(j,2)-spans(2)
         if (f.ge.0.0) go to 100
         f=-f/(spans(2)-spans(1))
         sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,1)
         go to 110
 100     f=f/(spans(3)-spans(2))
         sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,5)
 110  continue
      call smooth (n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
      return
      end
 