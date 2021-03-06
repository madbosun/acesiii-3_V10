      SUBROUTINE IDNTFY_SYMEQV_ATMS(PGRP, IORDGP, QREF, QTAR,
     &                              SYMEQVLNT) 
C
      IMPLICIT NONE 
C
      CHARACTER*4 PGRP 
      LOGICAL SYMEQVLNT, QSAME
      Integer IOrdGp, MAXORDER, I, J
      PARAMETER(MAXORDER = 8)
      Double precision QREF(3),QGEN(3, MAXORDER), QTAR(3)
C
C Given one set of Cartesian coordinates produces all symmetry
C equivalent coordinates and marks the symmetry unique ones. 
C
      DO 10 I = 1, IORDGP-1
         CALL SCOPY(3, QREF, 1, QGEN(1,I), 1)
10    CONTINUE
C
      IF (PGRP(1:3).EQ.'C2v')THEN
         QGEN(1,1) = -QGEN(1,1)
         QGEN(2,2) = -QGEN(2,2)
         QGEN(1,3) = -QGEN(1,3)
         QGEN(2,3) = -QGEN(2,3)
      ELSEIF(PGRP(1:3).EQ.'C2h')THEN
         QGEN(3,1) = -QGEN(3,1)
         QGEN(1,2) = -QGEN(1,2)
         QGEN(2,2) = -QGEN(2,2)
         QGEN(1,3) = -QGEN(1,3)
         QGEN(2,3) = -QGEN(2,3)
         QGEN(3,3) = -QGEN(3,3)
      ELSEIF(PGRP(1:3).EQ.'D2 ')THEN 
         QGEN(1,1) = -QGEN(1,1)
         QGEN(2,1) = -QGEN(2,1)
         QGEN(1,2) = -QGEN(1,2)
         QGEN(3,2) = -QGEN(3,2)
         QGEN(2,3) = -QGEN(2,3)
         QGEN(3,3) = -QGEN(3,3)
      ELSEIF(PGRP(1:3).EQ.'C2 ')THEN
         QGEN(1,1) = -QGEN(1,1)
         QGEN(2,1) = -QGEN(2,1)
      ELSEIF(PGRP(1:3).EQ.'C s')THEN
         QGEN(3,1) = -QGEN(3,1)
      ELSEIF(PGRP(1:3).EQ.'C i')THEN
         QGEN(1,1) = -QGEN(1,1)
         QGEN(2,1) = -QGEN(2,1)
         QGEN(3,1) = -QGEN(3,1)
      ELSEIF(PGRP(1:3).EQ.'D2h')THEN
C
C THIS COULD WELL BE WRONG!!!! CHECK!!!
C
         QGEN(1,1) = -QGEN(1,1)
         QGEN(2,2) = -QGEN(2,2)
         QGEN(1,3) = -QGEN(1,3)
         QGEN(2,3) = -QGEN(2,3)
         QGEN(3,4) = -QGEN(3,4)
         QGEN(1,5) = -QGEN(1,5)
         QGEN(3,5) = -QGEN(3,5)
         QGEN(2,6) = -QGEN(2,6)
         QGEN(3,6) = -QGEN(3,6)
         QGEN(1,7) = -QGEN(1,7)
         QGEN(2,7) = -QGEN(2,7)
         QGEN(3,7) = -QGEN(3,7)
      ENDIF
C
C Check and see whether the generated center is identical to the 
C target center. If that is the case then the reference center
C must be symmetry equivalent to the target center. 
C
      SYMEQVLNT = .FALSE. 
C
      DO 20 I = 1, IORDGP-1
C
       IF(QSAME(QTAR, QGEN(1,I))) SYMEQVLNT = .TRUE.
       IF (SYMEQVLNT) RETURN

20    CONTINUE
C 
      RETURN
      END
