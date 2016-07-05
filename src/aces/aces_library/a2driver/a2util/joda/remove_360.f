      SUBROUTINE REMOVE_360(DQ, TOTNOFBND, TOTNOFANG, TOTREDUNCO)
   
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER TOTNOFBND, TOTNOFANG, TOTREDUNCO
      DIMENSION DQ(TOTREDUNCO)
C
      TWOPI = (ATAN(DFLOAT(1))*DFLOAT(8))
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
C

      Print*, "The Two Phi, TOTNFBND, TOTNOFANG: ",
     &         TwOPI/DINVPI, TOTNOFBND, TOTNOFANG
      Write(6,*)
      Write(6,*) "Angles @Entry:"
      Write(6,*)
      Write(6,10) (DQ(I)/DINVPI, I = (TOTNOFBND + 1), TOTREDUNCO)
   10 Format (5(1X,F10.6)) 

C  
      DO IANG = TOTNOFANG+TOTNOFBND+1, TOTREDUNCO
C
            IF (DQ(IANG) .GE. TWOPI/2.0D0) THEN
                DQ(IANG) = DQ(IANG) - TWOPI
            ELSE IF (DQ(IANG) .LE. (-TWOPI/2.0D0)) THEN 
                DQ(IANG) = DQ(IANG) + TWOPI
            ENDIF
C
      ENDDO
C

      Print*, "The angles: @Exist"
      Write(6,*)
      Write(6,10) (DQ(I)/DINVPI, I = (TOTNOFBND + 1), TOTREDUNCO)

C
      RETURN
      END
