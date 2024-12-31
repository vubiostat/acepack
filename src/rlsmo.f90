! Refactored RLSMO and SMTH subroutines
MODULE SMOOTHING_MODULE
  IMPLICIT NONE
  CONTAINS

  ! Cross-validation for span selection.
  ! It evaluates a range of predefined spans selecting the optimal one
  ! using residual sums of squares
  SUBROUTINE RLSMO(x, y, w, span, dof, n, smo, rss, scrat)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: x(n), y(n), w(n)
    DOUBLE PRECISION, INTENT(INOUT) :: span
    DOUBLE PRECISION, INTENT(OUT) :: smo(n), rss, dof
    DOUBLE PRECISION, INTENT(OUT) :: scrat(n)
    DOUBLE PRECISION :: cvspan(6), cvrss(6), cvmin, penal, s0
    INTEGER :: k, idmin, i

    cvspan = (/0.3, 0.4, 0.5, 0.6, 0.7, 1.0/)
    penal = 0.01
    cvmin = 1.0d15
    idmin = 1

    ! Cross-validation for span selection
    IF (span == 0.0d0) THEN
      DO k = 1, 6
        CALL SMTH(x, y, w, cvspan(k), dof, n, 1, smo, s0, cvrss(k), scrat)
        IF (cvrss(k) < cvmin) THEN
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
    CALL SMTH(x, y, w, span, dof, n, 0, smo, s0, rss, scrat)
    DO i = 1, n
      smo(i) = smo(i) + s0
    END DO

  END SUBROUTINE RLSMO

  SUBROUTINE SMTH(x, y, w, span, dof, n, cross, smo, s0, rss, scrat)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, cross
    DOUBLE PRECISION, INTENT(IN) :: x(n), y(n), w(n), span
    DOUBLE PRECISION, INTENT(OUT) :: smo(n), s0, rss
    DOUBLE PRECISION, INTENT(OUT) :: scrat(n)
    DOUBLE PRECISION :: xbar, ybar, cov, var, sumw
    INTEGER :: is2, i

    ! Initialize variables
    s0 = 0.0d0
    rss = 0.0d0
    xbar = x(1)
    ybar = y(1)
    cov = 0.0d0
    var = 0.0d0
    sumw = w(1)

    IF (span < 1.0d0) THEN
      is2 = MAX(INT(span * n / 2), 1)

      DO i = 1, n
        ! Update neighborhood statistics
        CALL UPDATE_STATS(x, y, w, n, i, is2, xbar, ybar, cov, var, sumw)

        ! Calculate smoothed value
        IF (var > 0.0d0) THEN
          smo(i) = ybar + cov * (x(i) - xbar) / var
        ELSE
          smo(i) = ybar
        END IF
      END DO
    END IF

    ! Compute residual sum of squares
    DO i = 1, n
      rss = rss + (w(i) / sumw) * (y(i) - s0 - smo(i))**2
    END DO

  END SUBROUTINE SMTH

  SUBROUTINE UPDATE_STATS(x, y, w, n, i, is2, xbar, ybar, cov, var, sumw)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: x(n), y(n), w(n)
    DOUBLE PRECISION, INTENT(INOUT) :: xbar, ybar, cov, var, sumw
    INTEGER, INTENT(IN) :: n, i, is2

    INTEGER :: j
    DOUBLE PRECISION :: xin, yin, win

    ! Update statistics for the neighborhood
    DO j = MAX(1, i - is2), MIN(n, i + is2)
      xin = x(j)
      yin = y(j)
      win = w(j)
      xbar = (sumw * xbar + xin * win) / (sumw + win)
      ybar = (sumw * ybar + yin * win) / (sumw + win)
      cov = cov + win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw
      var = var + win * (xin - xbar)**2 * (sumw + win) / sumw
      sumw = sumw + win
    END DO

  END SUBROUTINE UPDATE_STATS

END MODULE SMOOTHING_MODULE
