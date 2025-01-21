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


! Cross-validation for span selection.
! It evaluates a range of predefined spans selecting the optimal one
! using residual sums of squares
! Cross-validation for span selection.
! It evaluates a range of predefined spans selecting the optimal one
! using residual sums of squares
SUBROUTINE rlsmo(x, y, w, span, dof, n, smo, rss, scratch)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: x(n)
  DOUBLE PRECISION, INTENT(IN)    :: w(n)
  DOUBLE PRECISION, INTENT(INOUT) :: y(n) ! FIXME: Why is this an INOUT ?
  DOUBLE PRECISION, INTENT(INOUT) :: span
  DOUBLE PRECISION, INTENT(OUT)   :: dof
  DOUBLE PRECISION, INTENT(OUT)   :: smo(n), rss
  DOUBLE PRECISION, INTENT(OUT)   :: scratch(n)
  
  DOUBLE PRECISION                :: cvspan(6), cvrss(6), cvmin, penal, s0
  INTEGER                         :: k, idmin

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
  smo(:) = smo(:) + s0

END SUBROUTINE rlsmo

