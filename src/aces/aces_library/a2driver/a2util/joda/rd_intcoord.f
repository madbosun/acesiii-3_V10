














































































































































































































































































































































































































































































































































      subroutine rd_intcoord(RINT)
      implicit none

      DOUBLE PRECISION DACOS, RINT
      INTRINSIC DACOS
      INTEGER LINBLNK
      EXTERNAL LINBLNK

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
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


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


C
      DOUBLE PRECISION dSmallNumber
      DIMENSION RINT(3*MXATMS)
      PARAMETER (dSmallNumber=5.361d-21)

      INTEGER        NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT
      COMMON /USINT/ NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT

      INTEGER        IFLAGS(100),IFLAGS2(500)
      COMMON /FLAGS/ IFLAGS,IFLAGS2

      DOUBLE PRECISION RTOA, dTmp
      INTEGER I, J, IBOMB, ISTP
      CHARACTER*(linelen) DUMSTR
      CHARACTER*5 CHTEST
      CHARACTER*1 czTab
      LOGICAL PRINTI

c ----------------------------------------------------------------------

C NEW (5/89) VARIABLE INPUT. NO LONGER HAVE TO INCLUDE REDUNDANCIES, AN
C UNIQUE COORDINATES CAN BE READ IN ANY ORDER. VARIABLE NAMES AND VAL
C ARE SEPARATED BY AN EQUALS (=) SIGN. MUCH BETTER THAN THE WAY IT WA
C DONE BEFORE. START BY FILLING R VECTOR WITH BOGUS VALUES.

      IF (IFLAGS(68).EQ.2) THEN

c Ajith Perera 12/2001.
C This block corresponds to the COORINATES=XYZ2INT option. To
C understand what this means, one needs to know the difference
C between COORDINATES=CARTESIAN and COORDINATES=XYZINT. In the first
C case, our input is the atomic symbol followed by the Cartesian
C coordinates of the each atom. There is no information about
C the internal coordinates. In the latter case, as one would
C have done for entering internal coordinates, one would enter
C the connectivity data in ZMAT format followed by (one blank
C in between) the Cartesian coordinates of the each atom in free
C format form.
C At present, this option does not understand the symmetry unless
C the user imposes it by ZMAT input. The advantage of having
C this capability is that one can use the Cartesian coordinates
C directly without first converting them to internal.
C In the future, we plan to extend this to include symmetry
C without user specfications. That way one can do optimiations
C with Cartesian input.
C
C Originaly, we had a routine called XTOR2.F which was pretty
C much doing what XTOR.F was doing except for the conversion of
C angles to radians. I obsoleted the XTOR2.F and extended the
C XTOR.F. The IPRT flag (don't know why it was named that)
C controls whether the conversion is needed (IPRT==1) or not
C (IPRT!=1).
         WRITE(LUOUT,*) '@RD_INTCOORD: Input coordinates are ',
     &      'transformed to internal representation.'
         READ(LUZ,*) (Q(I),I=1,NX)
         CALL XTOR(R,1)

      ELSE

         CALL GETREC(1,'JOBARC','FIRSTRUN',1,I)
         PRINTI = (I.NE.0)
         czTab  = achar(9)
         RTOA   = 180.D0/DACOS(-1.D0)
         IBOMB  = 0
         DO I = 1, NX
            R(I) = dSmallNumber
         END DO

         READ(LUZ,'(A)',END=987) DUMSTR
         do while (linblnk(dumstr(1:80)).ne.0)
            ISTP = INDEX(DUMSTR(1:),'=')
            I    = INDEX(DUMSTR(1:),'*')
            if ((i.ne.0).or.(istp.eq.0)) then
               write(*,*) '@RD_INTCOORD: Coordinate list is not ',
     &                    'terminated with a blank line.'
               call errex
            end if
            i = 1
            do while ((i.lt.istp).and.
     &                (dumstr(i:i).eq.' '.or.dumstr(i:i).eq.czTab))
               i = i + 1
            end do
            if (i.eq.istp) then
               write(*,*) '@RD_INTCOORD: Coordinate name is missing.'
               call errex
            end if
            j = i + 1
            do while ((j.lt.istp).and.
     &                (dumstr(j:j).ne.' '.and.dumstr(j:j).ne.czTab))
               j = j + 1
            end do
            CHTEST = DUMSTR(I:J-1)
            READ(DUMSTR(ISTP+1:),*) dTmp
            ISTP = 0
            I = 1
            do while (ISTP.eq.0.and.I.le.NXM6)
               IF (CHTEST.EQ.VARNAM(ISQUASH(I))) THEN
                  ISTP = 1
                  IF (RINT(I).NE.dSmallNumber) THEN
                     IBOMB = 1
                     WRITE(*,*) '@RD_INTCOORD: ',CHTEST,
     &                          ' multiply defined.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'TDA  ') THEN
                     RINT(I) = RTOA*DACOS(-1.0D0/3.0D0)
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable TDA set to exact tetrahedral angle.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'IHA  ') THEN
                     RINT(I) = RTOA*DACOS(1.D0/DSQRT(5.D0))
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable IHA set to exact icosahedral angle.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'OA   ') THEN
                     RINT(I) = RTOA*DACOS(-1.D0/3.D0)*0.5D0
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable OA set to 0.5 * exact tetrahedral angle.'
                  ELSE
                     R(I) = dTmp
                  END IF
               END IF
               I = I + 1
c           end do while (ISTP.eq.0.and.I.le.NXM6)
            end do
            IF (ISTP.EQ.0) THEN
               IBOMB = 1
               WRITE(*,*) '@RD_INTCOORD: ',CHTEST,' not in Z-matrix.'
            END IF
            READ(LUZ,'(A)',END=987) DUMSTR
c        end do while (linblnk(dumstr(1:80)).ne.0)
         end do
  987    CONTINUE

c      o copy unique coordinate definitions
         DO I = 1, NUNIQUE
            DO J = 1, NEQ(IUNIQUE(I))
               RINT(IEQUIV(I,J)) = RINT(IUNIQUE(I))
            END DO
         END DO

c      o check if all coordinates have been defined
         DO I = 1, NXM6
            IF (RINT(I).EQ.dSmallNumber) THEN
               IBOMB = 1
               WRITE(*,*)
     &               '@RD_INTCOORD: ',VARNAM(ISQUASH(I)),' is not ',
     &               'defined at ',ISQUASH(I)
            END IF
         END DO

c      o stop on all errors
         IF (IBOMB.EQ.1) CALL ERREX

      END IF

c   o unsquash R
      CALL USQUSH(R,NXM6)

      return
      end

