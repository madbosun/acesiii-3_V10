
      SUBROUTINE POLYPRT0(NATOMS,IINTFP,IPRINT,
     &                    LCOORD,LENERG,LGRAD,LHESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCOORD,LENERG,LGRAD,LHESS
      PARAMETER(MXATMS=250)
      DIMENSION  SCR(3*MXATMS + 3*MXATMS + 9*MXATMS*MXATMS)
      DIMENSION ISCR(  MXATMS)

      IF (NATOMS.GT.MXATMS) THEN
         WRITE(6,1000) NATOMS, MXATMS
 1000    FORMAT('  @POLYPRT0: Too many atoms. Given ',I3,' Max ',I3)
         CALL ERREX
      END IF

      I000 = 1
      I010 = I000 + 3*NATOMS
      I020 = I010 + 3*NATOMS
      I030 = I020 + 9*NATOMS*NATOMS
      CALL POLYPRT(SCR(I000),SCR(I010),SCR(I020),ISCR,
     &             NATOMS,IINTFP,IPRINT,
     &             LCOORD,LENERG,LGRAD,LHESS)

      RETURN
      END
