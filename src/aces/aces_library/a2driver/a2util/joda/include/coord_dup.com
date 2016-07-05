#ifndef _COORD_COM_
#define _COORD_COM_
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3, MXATMS), R(3, MAXREDUNCO/3), 
     &     NCON(3, MAXREDUNCO/3), NR(MXATMS),
     &     ISQUASH(MAXREDUNCO),IATNUM(MXATMS),ATMASS(MXATMS),
     &     IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end
#endif /* -COORD_COM_ */


