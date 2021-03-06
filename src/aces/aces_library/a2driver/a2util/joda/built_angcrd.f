      SUBROUTINE BULT_ANGCRD(CARTCOORD, BMATRX, ANGL, ICON1, ICON2, 
     &                       ICON3, IANGS, TOTREDNCO, NRATMS)
C         
C Setup the bond angle bending B-matrix elements. 
C Given r_i = r_ba and r_j = r_bc then CosG_ij = Sum (x_a-x_b)*(x_c-x_b)/R_i*R_j
C                                          x,y,z
C B(*,a,*,*)(x,y,z) = {CosG_ij*(x_a-x_b)/R_ba -(x_a-x_c)/R_ac}R_ba*SinG_ij
C B(*,*,*,c)(x,y,z) = {CosG_ij*(x_c-x_b)/R_bc -(x_a-x_b)/R_ab}R_bc*SinG_ij
C B(*,*,b,*)(x,y,z) = B(*,a,*,*) - B(*,*,*,c)
C
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION COSGAMA ,DISTBA ,DISTBC ,DISTAC, ANGL,
     &                 LEN_N2AC
      INTEGER TOTREDNCO
C
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO,3*NRATMS),
     &          VECBA(3), VECBC(3), VECBAD(3), VECN2AC(3),
     &          VECBN2A(3), VECBN2B(3)
C
      DATA ONE /1.0D0/
C
      COSGAMA = 0.0D00
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0  
C
      DISTBA = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
      DISTBC = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2))
      DISTAC = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON3 - 2))
C
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECBC, 1)
C 
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON1 - 2), 
     &         VECBA, 1)
C
      ANGL    = ANGLE(VECBC, VECBA, 3)*DINVPI
      SINGAMA = DSIN(ANGL)
      COSGAMA = DCOS(ANGL)
C
      IF (ANGL/DINVPI .LT. 175.0D0) THEN
C
         DO 20 IXYZ = 1, 3 
            BMATRX(IANGS, (3*ICON1 - 3) + IXYZ) = 
     &                    (COSGAMA*VECBA(IXYZ) - VECBC(IXYZ))
     &                                /(SINGAMA*DISTBA)
            BMATRX(IANGS, (3*ICON3 - 3) + IXYZ) = 
     &                   (COSGAMA*VECBC(IXYZ) - VECBA(IXYZ))
     &                                /(SINGAMA*DISTBC)
            BMATRX(IANGS, (3*ICON2 - 3) + IXYZ) = 
     &                - (BMATRX(IANGS, (3*ICON1 - 3) + IXYZ))
     &                - BMATRX(IANGS, (3*ICON3 - 3) + IXYZ)

 20      CONTINUE
C
      ELSE

          Write(6,*)
          Write(6, "(a)") "Entering the regular linear angle block"

C
          CALL NORMAL(VECBC, 3)
          CALL NORMAL(VECBA, 3)
          CALL CROSS(VECBC, VECBA, VECN2AC, 0)
C

          Write(6, *)
          Write(6, *) "The normal vector"
          Write(6, "(a, 3F10.4)") "VECN2AC", (VECN2AC(I), I=1,3)

C
          LEN_N2AC = DSQRT(DABS(XDOT(3, VECN2AC, 1, VECN2AC, 1)))
          IF (LEN_N2AC .LE. 1.0D-8) THEN
             VECN2AC(1) =  VECBC(2) + VECBC(3)
             VECN2AC(2) =  VECBC(1) - VECBC(3)
             VECN2AC(3) = -VECBC(1) - VECBC(2)
          ENDIF
C

          Write(6, *)
          Write(6, *) "The unit vectors"
          Write(6, "(a, 3F10.4)") "VECN2AC", (VECN2AC(I), I=1,3)

C
          LEN_N2AC = DSQRT(DABS(XDOT(3, VECN2AC, 1, VECN2AC, 1)))
          IF (LEN_N2AC .LE. 1.0D-8) THEN
             VECN2AC(1) =  VECBC(2) - VECBC(3)
             VECN2AC(2) =  VECBC(1) + VECBC(3)
             VECN2AC(3) =  VECBC(1) + VECBC(2)
          ENDIF
C

          Write(6, *)
          Write(6, *) "The unit vectors"
          Write(6, "(a, 3F10.4)") "VECN2AC", (VECN2AC(I), I=1,3)

C
          CALL NORMAL(VECN2AC, 3)
          CALL CROSS(VECBC, VECN2AC, VECBN2A, 1)
          CALL CROSS(VECN2AC, VECBA, VECBN2B, 1)
C

          Write(6, *) 
          Write(6, *) "The basis vctors for regular angles"
          Write(6, "(a, 3F10.4)") "VECBN2A", (VECBN2A(I), I=1,3)
          Write(6, "(a, 3F10.4)") "VECBN2B", (VECBN2B(I), I=1,3)
          Write(6,*)

C
          DO IXYZ = 1, 3
             BMATRX(IANGS, (3*ICON3 - 3) + IXYZ) = 
     &                                   VECBN2A(IXYZ)/DISTBC
             BMATRX(IANGS, (3*ICON2 - 3) + IXYZ) = 
     &                                  - VECBN2A(IXYZ)/DISTBC 
     &                                  - VECBN2B(IXYZ)/DISTBA
             BMATRX(IANGS, (3*ICON1 - 3) + IXYZ) = 
     &                                  VECBN2B(IXYZ)/DISTBA
          ENDDO
      ENDIF 
C
      RETURN
      END
      
