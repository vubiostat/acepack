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

! Original: https://stacks.stanford.edu/file/druid:gw754yg8889/ORION%20003.pdf
! J Friedman, W Stuetzle. Smoothing of Scatterplots, Stanford Project Orion, July 1982
SUBROUTINE SuperSmoother(x, y, w, span, dof, n, cross, smo, s0, rss, scratch)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)  :: x(n)   ! Ordered abscissa values
  DOUBLE PRECISION, INTENT(INOUT):: y(n)  ! Corresponding ordinate (response) values
  DOUBLE PRECISION, INTENT(IN)  :: w(n)   ! (optional) Weight for each (x,y) observation
  DOUBLE PRECISION, INTENT(IN)  :: span   ! Fractional span for residual smoothing
  DOUBLE PRECISION, INTENT(OUT) :: dof    ! Degrees of freedom
  INTEGER, INTENT(IN)           :: n      ! Number of observations (x,y)
  INTEGER, INTENT(IN)           :: cross  ! Cross validation
  DOUBLE PRECISION, INTENT(OUT) :: smo(n) ! Smoothed ordinate (response) variable
  DOUBLE PRECISION, INTENT(OUT) :: s0     ! Intercept
  DOUBLE PRECISION, INTENT(OUT) :: rss    ! Residual sum of squares
  DOUBLE PRECISION, INTENT(INOUT)  :: scratch(n) 
  
  INTEGER :: i, ibnew, ibold, is2, itnew, itold, j, jj, m0, ntie
  DOUBLE PRECISION :: xin, yin, win, xout, yout
  DOUBLE PRECISION :: wt, ispan
  DOUBLE PRECISION :: sumw, xbar, ybar, cov, var, r
  INTEGER :: fixeds, line

  line = 1
  fixeds = 1

  ! Check for invalid span value
  IF (span >= 1.0) THEN
    line = 0
  END IF

  ! Initialize the first values
  xbar = x(1)
  ybar = y(1)
  cov = 0.0
  var = 0.0
  sumw = w(1)

  IF (line == 1) THEN
    ! Main loop over the data
    DO i = 2, n
      CALL UpdateStats(x(i), y(i), w(i), xbar, ybar, cov, var, sumw, .FALSE.)
    END DO
  END IF

  DO i = 1, n
    IF (cross .EQ. 1) THEN
      CALL UpdateStats(x(i), y(i), w(i), xbar, ybar, cov, var, sumw, .TRUE.)
    END IF

    IF (var <= 0.0) THEN
      smo(i) = 0.0
    ELSE
      smo(i) = cov * (x(i) - xbar) / var
    END IF
  END DO

  s0 = ybar
  scratch(1) = cov / var
  itold = 1
  ibold = 1

  IF (cross == 0) THEN
    ! Initialize scratch
    DO i = 1, n
      scratch(i) = y(i)
    END DO

    ! Handle tied values
    DO i = 1, (n-1)
      m0 = i
      IF (x(i + 1) > x(i)) THEN
        EXIT
      END IF
    END DO

    IF (i == m0) THEN
      ntie = i - m0 + 1
      r = 0.0
      wt = 0.0

      DO jj = m0, i
        r = r + y(jj) * w(jj)
        wt = wt + w(jj)
      END DO

      r = r / wt
      DO jj = m0, i
        y(jj) = r
      END DO
    END IF
  END IF

  ! Final smoothing result
  ispan = n * span
  IF (fixeds /= 1) THEN
    is2 = INT(ispan / 2)
    IF (is2 < 1) is2 = 1
  END IF

  ! Final loop to update smoothing
  DO i = 1, n
    itnew = MIN(i + is2, n)
    ibnew = MAX(i - is2, 1)

    IF (itold >= itnew) THEN
      itold = itold + 1
      CALL UpdateStats(x(itold), y(itold), w(itold), xbar, ybar, cov, var, sumw, .FALSE.)
    END IF
  END DO

  ! Finalize
  DO i = 1, n
    y(i) = scratch(i)
  END DO

  ! Calculate the residual sum of squares (rss)
  rss = 0.0
  DO i = 1, n
    rss = rss + (w(i) / sumw) * (y(i) - s0 - smo(i))**2
  END DO

  ! Calculate the final degrees of freedom (dof)
  dof = dof + SUM(w) / sumw + SUM(w * (x - xbar)**2) / var

  RETURN
END SUBROUTINE SuperSmoother


SUBROUTINE UpdateStats(xin, yin, win, xbar, ybar, cov, var, sumw, is_remove)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)    :: xin, yin, win
  DOUBLE PRECISION, INTENT(INOUT) :: xbar, ybar, cov, var, sumw
  LOGICAL, INTENT(IN)             :: is_remove

  IF (is_remove) THEN
    ! If we are removing an element, subtract using (sumw - win)
    cov = cov - win * (xin - xbar) * (yin - ybar) * sumw / (sumw - win)
    var = var - win * (xin - xbar)**2 * sumw / (sumw - win)
    xbar = (sumw * xbar - win * xin) / (sumw - win)
    ybar = (sumw * ybar - win * yin) / (sumw - win)
    sumw = sumw - win
  ELSE
    ! Normal addition logic for updating stats
    xbar = (sumw * xbar + xin * win) / (sumw + win)
    ybar = (sumw * ybar + yin * win) / (sumw + win)
    cov = cov + win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw
    var = var + win * (xin - xbar)**2 * (sumw + win) / sumw
    sumw = sumw + win
  END IF
END SUBROUTINE UpdateStats
