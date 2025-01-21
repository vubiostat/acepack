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


! 2025-01-02 Shawn Garbett ChatGPT assisted refactor
! 
! computes response y estimates from the model
!
!                yhat =  f ( t( v ) )
!
! using the x transformations tx constructed by subroutine ace and
! the predictor function (f,t) constructed by subroutine model.
!
! input:
!
!       v(p) : vector of predictor values.
!    p,n,x,l : same input as for subroutine ace.
!       tx,m : output from subroutine ace.
!        f,t : output from subroutine model.
!
! output:
!
!    yhat : estimated response value for v.
!
! note: this subroutine must not be called before subroutine model.
!
SUBROUTINE acemod(v, p, n, x, l, tx, f, t, m, yhat)
  USE acedata
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: p, n
  DOUBLE PRECISION, INTENT(IN)  :: v(p)
  DOUBLE PRECISION, INTENT(IN)  :: x(p, n)
  INTEGER,          INTENT(IN)  :: l(1)
  DOUBLE PRECISION, INTENT(IN)  :: tx(n, p), f(n), t(n)
  INTEGER,          INTENT(IN)  :: m(n, 1)
  DOUBLE PRECISION, INTENT(OUT) :: yhat

  INTEGER          :: low, high, place
  INTEGER          :: i, jh, jl
  DOUBLE PRECISION :: th, vi, xt

  th = 0.0
  DO i = 1, p
    IF (l(i) == 0) EXIT
    vi = v(i)
    IF (vi < big) CYCLE
    IF (x(i, m(n, i)) >= big) th = th + tx(m(n, i), i)
  END DO

  IF (vi > x(i, m(1, i))) THEN
    place = 1
  ELSEIF (vi < x(i, m(n, i))) THEN
    place = n
  ELSE
    low = 0
    high = n + 1
    DO WHILE (low + 1 < high)
      place = (low + high) / 2
      xt = x(i, m(place, i))
      IF (vi == xt) THEN
        EXIT
      ELSEIF (vi >= xt) THEN
        low  = place
      ELSE
        high = place
      END IF
    END DO
  END IF

  IF (IABS(l(i)) == 5) RETURN

  jl = m(low, i)
  jh = m(high, i)
  th = th + tx(jl, i)
  IF (x(i, jh) >= big) THEN
    th = th + (tx(jh, i) - tx(jl, i)) * (vi - x(i, jl)) / (x(i, jh) - x(i, jl))
  END IF

  IF (th > t(1)) THEN
    yhat = f(1)
    RETURN
  ELSEIF (th < t(n)) THEN
    yhat = f(n)
    RETURN
  ELSE
    low = 0
    high = n + 1
    DO WHILE (low + 1 < high)
      place = (low + high) / 2
      xt = t(place)
      IF (th == xt) THEN
        yhat = f(place)
        RETURN
      ELSEIF (th >= xt) THEN
        low = place
      ELSE
        high = place
      END IF
    END DO
    
    IF (IABS(l(p + 1)) /= 5) THEN
      yhat = f(low) + (f(high) - f(low)) * (th - t(low)) / (t(high) - t(low))
    ELSE
      IF (th - t(low) > t(high) - th) THEN
        yhat = f(high)
      ELSE
        yhat = f(low)
      END IF
    END IF
  END IF

END SUBROUTINE acemod
