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
! supporting documentation. no representations are made about the
! suitability of this software for any purpose.  it is provided "as is"
! without express or implied warranty.
!______________________________________________________________________________

SUBROUTINE AddPoint(xin, yin, win, xbar, ybar, cov, var, sumw)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)    :: xin, yin, win
  DOUBLE PRECISION, INTENT(INOUT) :: xbar, ybar, cov, var, sumw

  xbar = (sumw * xbar + xin * win) / (sumw + win)
  ybar = (sumw * ybar + yin * win) / (sumw + win)
  cov = cov + win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw
  var = var + win * (xin - xbar)**2 * (sumw + win) / sumw
  sumw = sumw + win
END SUBROUTINE AddPoint

SUBROUTINE SubtractPoint(xin, yin, win, xbar, ybar, cov, var, sumw)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)    :: xin, yin, win
  DOUBLE PRECISION, INTENT(INOUT) :: xbar, ybar, cov, var, sumw

  cov = cov - win * (xin - xbar) * (yin - ybar) * sumw / (sumw - win)
  var = var - win * (xin - xbar)**2 * sumw / (sumw - win)
  xbar = (sumw * xbar - win * xin) / (sumw - win)
  ybar = (sumw * ybar - win * yin) / (sumw - win)
  sumw = sumw - win
END SUBROUTINE SubtractPoint

! Original: https://stacks.stanford.edu/file/druid:gw754yg8889/ORiOn%20003.pdf
! J Friedman, w Stuetzle. Smoothing of Scatterplots, Stanford Project Orion, July 1982
! This is the smoother used by ACE.
SUBROUTINE SuperSmoother(x,y,w,span,dof,n,cross,smo,s0,rss,scratch)
  IMPLICIT NONE
  INTEGER, INTENT(IN)           :: n      ! Number of observations (x,y)
  DOUBLE PRECISION, INTENT(IN)  :: x(n)   ! Ordered abscissa values
  DOUBLE PRECISION, INTENT(INOUT):: y(n)  ! Corresponding ordinate (response) values
  DOUBLE PRECISION, INTENT(IN)  :: w(n)   ! (optional) Weight for each (x,y) observation
  DOUBLE PRECISION, INTENT(IN)  :: span   ! Fractional span for residual smoothing
  DOUBLE PRECISION, INTENT(OUT) :: dof    ! Degrees of freedom
  INTEGER, INTENT(IN)           :: cross  ! Cross validation
  DOUBLE PRECISION, INTENT(OUT) :: smo(n) ! Smoothed ordinate (response) variable
  DOUBLE PRECISION, INTENT(OUT) :: s0     ! Intercept
  DOUBLE PRECISION, INTENT(OUT) :: rss    ! Residual sum of squares
  DOUBLE PRECISION, INTENT(INOUT)  :: scratch(n) 
  
  INTEGER :: i, ibnew, ibold, is2, itnew, itold, j, jj, m0, ntie
  DOUBLE PRECISION :: xin, yin, win, xout, yout
  DOUBLE PRECISION :: wt, ispan
  DOUBLE PRECISION :: sumw, xbar, ybar, cov, var, r
  INTEGER :: fixeds
  
  ! Initialize variables
  fixeds = 1
  xbar = x(1)
  ybar = y(1)
  cov = 0.0
  var = 0.0
  sumw = w(1)

  IF (span >= 1.0) THEN
    DO i=2, n
      CALL AddPoint(x(i), y(i), w(i), xbar, ybar, cov, var, sumw)
    END DO
    
    i = 1
    DO WHILE (i <= n)
      IF (cross == 1) CALL SubtractPoint(x(i), y(i), w(i), xbar, ybar, cov, var, sumw)
  
      IF (var <= 0.0) THEN
        smo(i) = 0.0
      ELSE
        smo(i) = cov * (x(i) - xbar) / var
      END IF
  
      IF (cross == 1) CALL AddPoint(x(i), y(i), w(i), xbar, ybar, cov, var, sumw)
      i = i + 1
    END DO
  
    s0 = ybar
    scratch(1) = cov / var
    dof = 1.0
    rss = sum( (w(:)/sumw) * (y(:) - s0 - smo(:)) ** 2 )
    RETURN
  END IF

  itold = 1
  ibold = 1
  dof = -1.0
  scratch(:) = y(:)

  IF (cross == 0) THEN 
    i = 0
    DO WHILE (i < n-1)
      i = i + 1
      m0 = i
      
      IF (x(i+1) <= x(i)) THEN
        DO 
          i = i + 1
          IF (i >= n) EXIT
        END DO
      END IF

      IF (i == m0) CYCLE
      
      ntie = i - m0 + 1
      
      r = 0.0
      wt = 0.0
      DO jj = m0, i
        j = jj
        r = r + y(j) * w(j)
        wt = wt + w(j)
      END DO

      r = r / wt
      DO j = m0, i
        y(j) = r
      END DO

    END DO
  END IF

  ispan = n * span
  IF (fixeds == 1) THEN
    is2 = int(ispan / 2)
    IF (is2 < 1) is2 = 1
  END IF

  DO i = 1, n
    itnew = min(i + is2, n)
    ibnew = max(i - is2, 1)
  
    DO WHILE (itold < itnew)
      itold = itold + 1
      CALL AddPoint(x(itold), y(itold), w(itold), xbar, ybar, cov, var, sumw)
    END DO

    DO WHILE (ibold > ibnew)
      ibold = ibold - 1
      CALL AddPoint(x(ibold), y(ibold), w(ibold), xbar, ybar, cov, var, sumw)
    END DO

    DO WHILE (itold > itnew)
      CALL SubtractPoint(x(itold), y(itold), w(itold), xbar, ybar, cov, var, sumw)
      itold = itold - 1
    END DO
  
    DO WHILE (ibold < ibnew)
      CALL SubtractPoint(x(ibold), y(ibold), w(ibold), xbar, ybar, cov, var, sumw)
      ibold = ibold + 1
    END DO
    
    IF (cross == 1) CALL SubtractPoint(x(i), y(i), w(i), xbar, ybar, cov, var, sumw)

    IF (var <= 0) THEN
      smo(i) = ybar
      dof = dof + w(i) / sumw
    ELSE 
      smo(i) = ybar + cov * (x(i) - xbar) / var
      dof = dof + w(i) / sumw + (w(i) * (x(i) - xbar) ** 2) / var
    END IF

    IF (cross == 1) CALL AddPoint(x(i), y(i), w(i), xbar, ybar, cov, var, sumw)

  END DO

  y(:) = scratch(:)

  IF (cross == 0) THEN
    i = 0
    
    DO WHILE (i < n-1)
      i  = i + 1
      m0 = i  
    
      DO WHILE (x(i+1) <= x(i) .and. i < n-1)
        i = i + 1
      END DO

      IF (i /= m0) THEN
        ntie = i - m0 + 1
        wt = SUM(w(m0:i))
        r = SUM(smo(m0:i) * w(m0:i))/wt
        smo(m0:i) = r
      END IF
    END DO
  END IF
  
  ybar   = sum(w(:)*y(:))
  sumw   = sum(w(:))
  ybar   = ybar / sumw
  smo(:) = smo(:) - ybar
  s0     = ybar
  rss    = sum( (w(:)/sumw) * (y(:) - s0 - smo(:)) ** 2 )
      
END SUBROUTINE SuperSmoother
