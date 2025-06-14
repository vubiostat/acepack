  !----------------------------------------------------------------------------
 !
! This file is part of acepack.
!
! Copyright 1985,2007 Jerome H. Friedman
! Copyright 2016,2025 Shawn Garbett, Vanderbilt University
!
! Permission to use, copy, modify, distribute, and sell this software and
! its documentation for any purpose is hereby granted without fee,
! provided that the above copyright notice appear in all copies and that
! both that copyright notice and this permission notice appear in
! supporting documentation. No representations are made about the
! suitability of this software for any purpose.  It is provided "as is"
! without express or implied warranty.
!______________________________________________________________________________


! Computes response predictive  function f for the model yhat = f(t),
! where
!                                        p
!            f(t) = e(y : t),     t =   sum  tx<i> ( x<i> )
!                                       i=1
! using the x transformations tx constructed by subroutine ace.
! if y is a categorical variable (classification) then
!                                -1
!                       f(t) = ty  (t).
! input:
!
!    p,n,y,w,l : same input as for subroutine ace.
!    tx,ty,m,z : output from subroutine ace.
!
! output:
!
!    f(n),t(n) : input for subroutine acemod.
!
! note: this subroutine must be called before subroutine acemod.
!

SUBROUTINE model(p, n, y, w, l, tx, ty, f, t, m, z)
  USE acedata
  IMPLICIT NONE

  ! Input/Output arguments
  INTEGER, INTENT(IN)           :: p, n
  INTEGER, INTENT(INOUT)        :: l(1), m(n, 1)
  DOUBLE PRECISION, INTENT(IN)  :: y(n), w(n), tx(n, p), ty(n)
  DOUBLE PRECISION, INTENT(OUT) :: f(n), t(n), z(n, 12)

  ! Local variables
  INTEGER          :: j, k, j1, j2, pp1

  pp1 = p + 1

  IF (ABS(l(pp1)) == 5) THEN
    t(:) = ty(:)
    m(:, pp1) = reshape([(j, j=1,n)], shape=[n]) ! 1:n
  ELSE
    DO j = 1, n
      t(j) = sum(tx(j, :))
      m(j, pp1) = j
    END DO
  END IF

  CALL sort(t, m(:, pp1), 1, n)

  ! Loop for populating z
  DO j = 1, n
    k = m(j, pp1)
    z(j, 2) = w(k)
    IF (y(k) >= big) CYCLE ! Skip updating z(j, 1)
    z(j, 1) = y(k)
  END DO
  
  j1=j
  j2=j1
  
  DO WHILE (y(m(j1,pp1)) >= big)
    j1 = j1 - 1
    IF (j1 < 1) EXIT
  END DO
  
  DO WHILE (y(m(j2,pp1)) >= big)
    j2 = j2 + 1
    IF (j2 > n) EXIT
  END DO

  IF (j1 < 1) THEN
    k=j2
  ELSE IF (j2 > n) THEN
    k=j1
  ELSE IF (t(j)-t(j1).ge.t(j2)-t(j)) THEN
    k=j2
  ELSE
    k=j1
  END IF

  z(j, 1) = y(m(k, pp1))
  t(j) = t(k)

  ! Smoothness check
  IF (ABS(l(pp1)) == 5) THEN
    f(:) = z(:, 1)
  ELSE
    CALL smothr(1, n, t, z, z(:, 2), f, z(:, 6))
  END IF

  RETURN
END SUBROUTINE model
