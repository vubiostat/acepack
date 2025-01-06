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


! Puts into a the permutation vector which sorts v into
! increasing order. Only elements from ii to jj are considered.
! arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements.
!     
! This is a modification of CACM algorithm #347 which is a modified
! Hoare's quicksort.
!
! Richard C. Singleton. (1969). Algorithm 347: an efficient algorithm for
! sorting with minimal storage [M1]. _Communications of the ACM_, 23(3), 185-7.
! doi: 10.1145/362875.362901.
  
SUBROUTINE sort (v,a,ii,jj)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)    :: ii, jj
  DOUBLE PRECISION, INTENT(INOUT) :: v(*)
  INTEGER,          INTENT(INOUT) :: a(jj)

  INTEGER           iu(20), il(20)
  INTEGER           t,tt,ij,j,k,l
  INTEGER           m,i
  DOUBLE PRECISION  vt,vtt
  LOGICAL           continue
      
  m=1
  i=ii
  j=jj
  
  DO 
    DO 
      IF (i < j) THEN 
        k=i
   
        ij=(j+i)/2
        t=a(ij)
        vt=v(ij)
    
        IF (v(i) > vt) THEN
          a(ij)=a(i)
          a(i)=t
          t=a(ij)
          v(ij)=v(i)
          v(i)=vt
          vt=v(ij)
        END IF
    
        l=j
        
        IF (v(j) < vt) THEN
          a(ij)=a(j)
          a(j)=t
          t=a(ij)
          v(ij)=v(j)
          v(j)=vt
          vt=v(ij)
          
          IF (v(i) > vt) THEN
            a(ij)=a(i)
            a(i)=t
            t=a(ij)
            v(ij)=v(i)
            v(i)=vt
            vt=v(ij)
          END IF
        END IF
        
        DO
          DO 
            l=l-1
            IF (v(l) <= vt) EXIT
          END DO
          
          tt=a(l)
          vtt=v(l)
          
          DO
            k=k+1
            IF (v(k) >= vt) EXIT
          END DO
          
          IF (k > l) EXIT
          
          a(l)=a(k)
          a(k)=tt
          v(l)=v(k)
          v(k)=vtt
          
        END DO
        
        IF (l-i > j-k) THEN
          il(m)=i
          iu(m)=l
          i=k
          m=m+1
          
          if (j-i > 10) CYCLE
          IF (i /= ii) EXIT
          CYCLE
        END IF
        
        il(m)=k
        iu(m)=j
        j=l
        m=m+1
  
        if (j-i > 10) CYCLE
        IF (i /= ii) EXIT
        CYCLE
        
      END IF
        
      m=m-1
      IF (m.eq.0) RETURN
      i=il(m)
      j=iu(m)
        
      IF (j-i > 10) CYCLE
      IF (i /= ii) EXIT
    END DO
      
    i=i-1
    
    continue = .FALSE.
    DO
      i=i+1
      IF (i == j) THEN
        continue = .TRUE.
        EXIT
      END IF
  
      t=a(i+1)
      vt=v(i+1)
      if (v(i).le.vt) CYCLE
      k=i
      DO
        a(k+1)=a(k)
        v(k+1)=v(k)
        k=k-1
        IF (vt >= v(k)) EXIT
      END DO
      
      a(k+1)=t
      v(k+1)=vt
    END DO
    
    IF (.not. continue) EXIT
    
  END DO
        
END SUBROUTINE sort