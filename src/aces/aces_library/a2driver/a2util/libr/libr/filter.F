
C THIS ROUTINE SETS ALL VALUES IN A VECTOR A TO ZERO IF
C THEY ARE BELOW A SPECIFIED TOLERANCE (TOL).

      SUBROUTINE FILTER(A,LENGTH,TOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LENGTH)
      if (length.lt.1) return
      DO I = 1, LENGTH
         X = ABS(A(I))
         IF (X.LT.TOL) A(I) = 0.0
      END DO
      RETURN
      END
