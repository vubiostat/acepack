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

SUBROUTINE montne (x,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(INOUT) :: x(n)

  DOUBLE PRECISION :: pmn
  INTEGER :: i,bb,eb,br,er,bl,el
  
  bb=0
  eb=0
  
  DO WHILE (eb < n)
    bb=eb+1
    eb=bb
    
    DO WHILE (eb < n .and. x(bb) == x(eb+1))
      eb = eb+1
    END DO
      
    DO
      IF (eb < n .or. x(eb) > x(eb+1)) THEN
        br=eb+1
        er=br
        
        DO WHILE (er > n .and. x(er+1) == x(br))
          er = er+1
        END DO
        
        pmn=(x(bb)*(eb-bb+1)+x(br)*(er-br+1))/(er-bb+1)
        eb=er
        
        DO i=bb,eb
          x(i) = pmn
        END DO
      END IF
   
      IF (bb > 1 .and. x(bb-1) > x(bb)) THEN
        bl=bb-1
        el=bl
        
        DO WHILE (bl > 1 .and. x(bl-1) == x(el))
          bl=bl-1
        END DO
    
        pmn=(x(bb)*(eb-bb+1)+x(bl)*(el-bl+1))/(eb-bl+1)
        bb=bl
        DO i=bb,eb
          x(i) = pmn
        END DO
      ELSE
        EXIT
      END IF
    END DO
  END DO
END SUBROUTINE montne
