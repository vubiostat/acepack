! Cross-validation for span selection.
! It evaluates a range of predefined spans selecting the optimal one
! using residual sums of squares
! Cross-validation for span selection.
! It evaluates a range of predefined spans selecting the optimal one
! using residual sums of squares
SUBROUTINE RLSMO(x, y, w, span, dof, n, smo, rss, scratch)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)    :: x(n)
  DOUBLE PRECISION, INTENT(IN)    :: w(n)
  DOUBLE PRECISION, INTENT(INOUT) :: y(n) ! FIXME: Why is this an INOUT ?
  DOUBLE PRECISION, INTENT(INOUT) :: span
  DOUBLE PRECISION, INTENT(OUT)   :: dof
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), rss
  DOUBLE PRECISION, INTENT(OUT)   :: scratch(n)
  
  DOUBLE PRECISION                :: cvspan(6), cvrss(6), cvmin, penal, s0
  INTEGER                         :: k, idmin, i

  cvspan = (/0.3, 0.4, 0.5, 0.6, 0.7, 1.0/)
  penal = 0.01
  cvmin = 1.0d15
  idmin = 1

  ! Cross-validation for span selection
  IF (span == 0.0d0) THEN
    DO k = 1, 6
      CALL SuperSmoother(x, y, w, cvspan(k), dof, n, 1, smo, s0, cvrss(k), scratch)
      IF (cvrss(k) <= cvmin) THEN
        cvmin = cvrss(k)
        idmin = k
      END IF
    END DO
    
    span = cvspan(idmin)

    IF (penal > 0.0d0) THEN
      cvmin = (1.0d0 + penal) * cvmin
      DO k = 6, 1, -1
        IF (cvrss(k) <= cvmin) THEN
          span = cvspan(k)
          EXIT
        END IF
      END DO
    END IF 
  END IF

  ! Final smoothing
  CALL SuperSmoother(x, y, w, span, dof, n, 0, smo, s0, rss, scratch)
  DO i = 1, n
    smo(i) = smo(i) + s0
  END DO

END SUBROUTINE RLSMO

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

! Original: https://stacks.stanford.edu/file/druid:gw754yg8889/ORION%20003.pdf
! J Friedman, W Stuetzle. Smoothing of Scatterplots, Stanford Project Orion, July 1982
SUBROUTINE SuperSmoother(x, y, w, span, dof, n, cross, smo, s0, rss, scratch)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)  :: x(n)   ! Ordered abscissa values
  DOUBLE PRECISION, INTENT(INOUT):: y(n)   ! Corresponding ordinate (response) values
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
END SUBROUTINE
