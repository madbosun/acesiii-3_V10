      SUBROUTINE OCCUPY(NIRREP,NBASIR,NBASTOT,EVAL,ISCR,NOCC,ISPIN)
C
C DETERMINES THE OCCUPANCY FOR AN SCF DETERMINANT BASED ON
C  AN EXISTING SET OF EIGENVALUES
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EVAL(NBASTOT),ISCR(NBASTOT,2),NSPIN(2)
      DIMENSION NBASIR(NIRREP), Nocc(8,2)
C
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
C
      IONE=1
      CALL GETREC(20,'JOBARC','NMPROTON',IONE,NPROTON)
      ICHARG=IFLAGS(28)
      IMULT=IFLAGS(29)
C
C DETERMINE ALPHA AND BETA OCCUPANCIES
C
      NUMEL=NPROTON-ICHARG
      IALPEX=IMULT-1
      NRHS=NUMEL-IALPEX
      NSPIN(2)=NRHS/2
      NSPIN(1)=NSPIN(2)+IALPEX
C
      IBASIS=0
      DO 10 I=1,NIRREP
       DO 20 J=1,NBASIR(I)
         IBASIS=IBASIS+1
         ISCR(IBASIS,1)=IBASIS
         ISCR(IBASIS,2)=I
20      CONTINUE
10     CONTINUE
C
       CALL IZERO(NOCC(1,ISPIN),8)
       CALL PIKSR2(NBASTOT,EVAL,ISCR)
       DO 40 IOCC=1,NSPIN(ISPIN)
        IPOS=ISCR(IOCC,1)
        IRREP=ISCR(IPOS,2)
        NOCC(IRREP,ISPIN)=NOCC(IRREP,ISPIN)+1
40     CONTINUE
30    CONTINUE
C
      RETURN
      END
