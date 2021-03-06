      SUBROUTINE SYMMET4(A,NDIM,NUM)
C
C   THIS ROUTINE SYMMETRIZES A GIVEN MATRIX A
C
C     A(PQ) = 1/2 ( A(PQ)+A(QP))
C
C   WHERE A IS A SYMMETRY PACKED MATRIX AND
C   NUM THE CORRESPONDING POPULATION VECTOR
C   NDIM IS THE LEADING DIMENSION OF THE MATRIX A
C
C   THE SYMMETRIZATION IS HERE COMPLETELY DONE IN PLACE
C
C   THIS ROUTINE IS REQUIRED IN THE OUT-CORE-SYMMETRIZATION
C   ROUTINE AND CALLED FROM SYMMET3
C
CEND
C
C CODED AUGUST/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NDIM,NUM)
      DATA HALF /0.5D+0/
C
C SYMMETRIZE FISRT HALF OF A
C
       DO 100 I=2,NUM
CDIR$ IVDEP
*VOCL LOOP,NOVREC
       DO 100 J=1,I-1
        A(I,J)=HALF*(A(I,J)+A(J,I))
100    CONTINUE
C
C FILL SECOND HALF OF A
C
       DO 200 I=2,NUM
CDIR$ IVDEP
*VOCL LOOP,NOVREC
       DO 200 J=1,I-1
        A(J,I)=A(I,J)
200    CONTINUE
C
      RETURN
      END
