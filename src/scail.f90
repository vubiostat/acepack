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


! 2025-01-02 Shawn Garbett
!   Refactored with ChatGPT assistance. 
!
! 2016-10-07 Shawn Garbett
!   Refactor to insure initialized variable h, and no division by zero.
!
! Note: The original function was named "scale", but this is now part of the
!       Fortran 95 namespace. So this was changed to "scail"
!
SUBROUTINE scail(p, n, w, sw, ty, tx, eps, maxit, r, sc)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: p, n, maxit
  DOUBLE PRECISION, INTENT(IN)    :: w(n), ty(n), sw, eps
  DOUBLE PRECISION, INTENT(INOUT) :: tx(n, p)
  DOUBLE PRECISION, INTENT(OUT)   :: r(n), sc(p, 5)

  INTEGER          :: i, j, iter, nit
  DOUBLE PRECISION :: s, h, t, u, gamma, delta, v
  DOUBLE PRECISION :: residual, product

  ! Initialization
  sc(:,1) = 0.0
  nit = 0

  DO
    nit = nit + 1   ! nit is iteration number
    
    ! Store the previous value of sc(i, 1) for convergence check
    sc(:,5) = sc(:,1)

    h = 1.0 ! Initialize h to avoid uninitialized variable warning

    ! Iterate over all steps for the current iteration
    DO iter = 1, p

      ! Compute residuals: r(j) = (ty(j) - sum(sc(:,1) * tx(j,:))) * w(j)
      DO j = 1, n
        residual = ty(j) - DOT_PRODUCT(sc(:,1), tx(j,:))
        r(j) = residual * w(j)
      END DO

      ! Update sc(:,2) with the computed values
      sc(:,2) = -2.0 * (MATMUL(r, tx) / sw)

      ! Calculate s for next iteration and avoid zero denominator
      s = SUM(sc(:,2)**2)

      ! Ensure h gets initialized with s
      IF (iter == 1 .OR. h <= 0.0) THEN
        h = s
      END IF

      ! If sum of sc(i,2)^2 is zero, exit the loop
      IF (s <= 0.0) EXIT

      ! Update sc(:,3) based on gamma and previous values
      gamma = s / h
      sc(:,3) = -sc(:,2) + gamma * sc(:,4)

      h = s

      ! Calculate new delta and update sc(i,1) and sc(i,4)
      s = 0.0
      t = 0.0
      DO j = 1, n
        product = DOT_PRODUCT(sc(:,3), tx(j,:))
        s = s + product * r(j)
        t = t + w(j) * product**2
      END DO
      delta = s / t

      ! Update sc(i,1) and sc(i,4) based on delta
      sc(:,1) = sc(:,1) + delta * sc(:,3)
      sc(:,4) = sc(:,3)
    END DO

    ! Check for convergence
    v = MAXVAL(ABS(sc(:,1) - sc(:,5)))
    IF (v < eps .OR. nit >= maxit) EXIT
  END DO

  ! Update tx(j,i) based on the final sc(i,1)
  DO i = 1, p
    tx(:,i) = sc(i,1) * tx(:,i)
  END DO

  RETURN
END SUBROUTINE scail
