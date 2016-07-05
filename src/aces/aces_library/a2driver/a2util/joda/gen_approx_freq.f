













































































































































































































































































































































































































































































































































      SUBROUTINE GEN_APPROX_FREQ(SCRATCH, BMATRIX, TMP, HES_INTACT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)


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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end









c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C
      LOGICAL XYZIN, NWFINDIF
      CHARACTER*1 FLAG(3*MXATMS)
      INTEGER TOTREDNCO
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DOUBLE PRECISION BMATRIX(NXM6*3*NATOMS), TMP(NXM6*3*NATOMS),
     &                 SCRATCH(9*NATOMS*NATOMS),
     &                 HES_INTACT(3*NATOMS, 3*NATOMS)

C Transform the internal coordinate Hessian to the Cartesian Hessian in
C order to carry out mass weighing (it is nearly impossible to mass
C weigh internal Hessians). Notice the assumption made in the
C transformation: the gradient is zero, hence no derivative of the B
C matrix is required. Since this is done only at the stationary points,
C the assumption (or the condition) is satisfied and the transformation
C is valid. Ajith Perera 05/2005.

CSSS      CALL OUTPUT(HES_INTACT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
      IF (Iflags2(5) .ge. 3) THEN
C
C It is more than safer to always use 3*NATOMS explicitly when
C redundent internal optimizations are performed. Ajith perera,
C 01/2007.

         CALL GETREC(20, 'JOBARC','REDNCORD', 1, TOTREDNCO)
         LENGTH_BMAT = 3*NATOMS*TOTREDNCO
         CALL GETREC(-1, 'JOBARC', 'BMATRIXT', LENGTH_BMAT*IINTFP,
     &              TMP)
         CALL ZERO(SCRATCH, 9*NATOMS*NATOMS)
         CALL XGEMM("N", "N", 3*NATOMS, TOTREDNCO, TOTREDNCO, 1.0D0,
     &              TMP, 3*NATOMS, HES_INTACT, TOTREDNCO, 0.0D0,
     &              SCRATCH, 3*NATOMS)
         CALL TRANSP(TMP, BMATRIX, TOTREDNCO, 3*NATOMS)
         CALL ZERO(HES_INTACT, 9*NATOMS*NATOMS)
         CALL XGEMM("N", "N", 3*NATOMS, 3*NATOMS, TOTREDNCO, 1.0D0, 
     &               SCRATCH, 3*NATOMS, BMATRIX, TOTREDNCO, 0.0D0, 
     &               HES_INTACT, 3*NATOMS)
C
      ELSE
C
         CALL GETREC(20, 'JOBARC', 'BMATRIX ', IINTFP*NXM6*NX,
     &               BMATRIX)
CSSS         CALL OUTPUT(BMATRIX, 1, NXM6, 1, NX, NXM6, NX, 1)
C
         CALL TRANSP(BMATRIX, TMP, NX, NXM6)
         CALL ZERO(SCRATCH, NX*NX)
         CALL XGEMM("N", "N", NX, NXM6, NXM6, 1.0D0, TMP, NX,
     &               HES_INTACT, NXM6, 0.0D0, SCRATCH, NX)
         CALL ZERO(HES_INTACT, NX*NX)
         CALL XGEMM("N", "N", NX, NX, NXM6, 1.0D0, SCRATCH, NX,
     &               BMATRIX, NXM6, 0.0D0, HES_INTACT, NX)
C
CSSS         CALL OUTPUT(HES_INTACT, 1, NX, 1, NX, NX, NX, 1)
      ENDIF
C
C Mass weight the Cartesian Hessian
C
      DO I = 1, 3*NATOMS
         DO J = 1, 3*NATOMS
            FACT = DSQRT(ATMASS(1+(I-1)/3)*ATMASS(1+(J-1)/3))
            IF (FACT .LT. 1.0D-3)THEN
                HES_INTACT(I,J) = 0.D0
            ELSE
                HES_INTACT(I,J) = HES_INTACT(I,J)/FACT
            ENDIF
         ENDDO
      ENDDO
C
C Diagonalize the mass-weighted Hessian and tag the negative
C eigenvalues with "i" to indicate an imaginary frequency and
C convert them to the correspnding frequencies in cm-1 units.
C
      CALL EIG(HES_INTACT, SCRATCH, 3*NATOMS, 3*NATOMS, 0)
C
      DO I = 1, 3*NATOMS
         FLAG(I) = " "
         IF (HES_INTACT(I,I) .LT. 0.0D0) FLAG(I) = "i"
         HES_INTACT(I,I) = DSQRT(DABS(HES_INTACT(I,I)))*5.14048D03
         SCRATCH(I) = HES_INTACT(I,I)
      END DO
C
C Print the frequencies along with a warning message.
C
      WRITE(LUOUT,10)
   10 FORMAT(3X, "Frequencies of the updated Hessian at convergence:")
      WRITE(LUOUT,*)
C
      WRITE(LUOUT,'(11(F6.0,A))')
     &     (SCRATCH(IMODE), FLAG(IMODE), IMODE=1,3*NATOMS)
      WRITE(LUOUT,*)
      WRITE(LUOUT,20)
   20 FORMAT(3X, "Warning: Frequencies based on the updated Hessian",
     &" are not very accurate.",/,12X, "A separate frequency",
     &" calculation is required to get correct values.")
      WRITE(LUOUT,*)
C
      RETURN
      END

