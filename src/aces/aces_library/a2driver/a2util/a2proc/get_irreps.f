C
C DRIVER FOR DETERMINATION OF ORBITAL SYMMETRIES. From ACES II. 
C 
      SUBROUTINE GET_IRREPS(SCFVEC,SCFEVAL,ICORE,MAXCOR,NBAS,NBASX,
     &                      ISPIN, NOCC, IUHF)
C
      IMPLICIT INTEGER (A-Z)

      DIMENSION ICORE(MAXCOR), IREPS(9), NBFIRR(8), NOCC(8,2)
      DOuble Precision Scfvec(Nbasx*Nbas), Scfeval(Nbas)
      CHARACTER*1 SP(2)
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



C
      DATA IONE /1/
      DATA SP /'A','B'/

C FOR SPHERICAL HARMONICS, DIMENSION OF COMPUTATIONAL BASIS
C IS SMALLER THAN THE CARTESIAN SET.  READ THAT VALUE NOW.

      CALL GETREC(-1,'JOBARC','NATOMS  ',IONE,NATOMS)
      CALL GETREC(-1,'JOBARC','COMPNIRR',IONE,NIRREP)
      CALL GETREC(-1,'JOBARC','FULLORDR',IONE,IORDGP)
      CALL GETREC(-1,'JOBARC','NUMBASIR',NIRREP, NBFIRR)

      CALL IZERO(ICORE(1),MAXCOR)
      NSIZE = MAX(NBAS,NATOMS,NIRREP,NBASX)

      I020 = 1
      I030 = I020 + IINTFP*MAX(NBASX*10,NATOMS)
      I040 = I030 + MAX(NBASX+MOD(NBASX,2),IINTFP*3*NATOMS)
      I060 = I040 + MAX(NBASX+MOD(NBASX,2),NATOMS+MOD(NATOMS,2))
      I080 = I060 + MAX(1080,NATOMS+MOD(NATOMS,2),IINTFP*NBASX*5)
      I090 = I080 + MAX(NATOMS*5+MOD(NATOMS,2),IINTFP*NBASX)
      I100 = I090 + NATOMS*5+MOD(NATOMS,2)
      I070 = I100 + IINTFP*IORDGP
      I110 = I070 + MAX(2000,NATOMS+MOD(NATOMS,2),NBASX+MOD(NBASX,2),
     &                  NBASX*NBASX*IINTFP)

      CALL IRRORB(SCFVEC,SCFEVAL,ICORE(I020),ICORE(I030),
     &            ICORE(I040),NBAS,                   ICORE(I060),
     &            ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     &            NATOMS,     ISPIN,      NBASX)
C
C Write a Table of Eigenvalues and their symmetry. 
C
      IREPS(1) = 1
      DO IRREP = 1, NIRREP
         IREPS(IRREP+1) = IREPS(IRREP) + NBFIRR(IRREP)
      ENDDO 

      I020 = 1
      I030 = I020 + NBASX
      I040 = I030 + NBASX
  
      CALL EVLOUT(SCFEVAL, ICORE(I020), ICORE(I030), IREPS, NBAS,
     &            NIRREP, NOCC, ISPIN)

      RETURN
      END

