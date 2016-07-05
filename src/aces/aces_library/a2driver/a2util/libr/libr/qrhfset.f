      SUBROUTINE QRHFSET
C
C  THIS ROUTINE SETS UP THE QRHFINF COMMON BLOCK, WHICH IS
C  VERY USEFUL FOR THE EVALUATION OF THE QRHF Z-VECTOR EQUATIONS
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER POP,VRT,POPRHF,VRTRHF,POPDOC,VRTDOC,DIRPRD
C



c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN

      PARAMETER (MAXBASFN=1000)



c maxbasfn.par : end


      DIMENSION ISET(255),SCR(2*MAXBASFN)
C
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2) 
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,N2A,
     &               NUMISCF,NUMASCF
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
C
      IONE=1
C
      CALL GETREC(20,'JOBARC','QRHFTOT ',IONE,NMOD)
      if (nmod.gt.255) then
         print *, '@QRHFSET: Assertion failed.'
         print *, '          iset dimension = 255'
         print *, '          nmod = ',nmod
         call errex
      end if
      CALL GETREC(20,'JOBARC','QRHFIRR ',NMOD,ISET)
c
c this is a temporary kludge until qrhfx, etc routines are
c available
c
c      if(iset(1).lt.0)iflags(32)=-iset(1)
c      if(iset(1).gt.0)iflags(33)=iset(1)
c Uncommenting the next line will break QRHF first-order props.
c They were tested numerically just to make sure that they are correct.
c 04/2006, Ajith Perera.
c      iflags(38)=0
C
      CALL ICOPY(NIRREP,POP(1,1),1,POPRHF,1)
      CALL ICOPY(NIRREP,VRT(1,2),1,VRTRHF,1)
C 
      DO 10 I=1,NMOD
       INDEX=ISET(I)
       IF(INDEX.GT.0)THEN
        POPRHF(INDEX)=POPRHF(INDEX)-1
       ELSE
        VRTRHF(-INDEX)=VRTRHF(-INDEX)-1
       ENDIF
10    CONTINUE
C
      NAI=0
      N1I=0
      N2A=0
      NUMISCF=0
      NUMASCF=0
      DO 20 IRREP=1,NIRREP
       POPDOC(IRREP)=POP(IRREP,2)
       VRTDOC(IRREP)=VRT(IRREP,1)
       NOSH1 (IRREP)=POPRHF(IRREP)-POP(IRREP,2)
       NOSH2 (IRREP)=POP(IRREP,1)-POPRHF(IRREP)
       NAI=NAI+POPRHF(IRREP)*VRTRHF(IRREP)
       N1I=N1I+NOSH1(IRREP)*POPDOC(IRREP)
       N2A=N2A+NOSH2(IRREP)*VRTDOC(IRREP)
       NUMISCF=NUMISCF+POPRHF(IRREP)
       NUMASCF=NUMASCF+VRTRHF(IRREP)
20    CONTINUE
      NBAS=NUMISCF+NUMASCF
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@QRHFSET-I, Contents of QRHFINF common block:')
       WRITE(6,1001)(POPRHF(I),I=1,NIRREP)
       WRITE(6,1001)(VRTRHF(I),I=1,NIRREP)
       WRITE(6,1001)(NOSH1(I),I=1,NIRREP)
       WRITE(6,1001)(NOSH2(I),I=1,NIRREP)
       WRITE(6,1001)(POPDOC(I),I=1,NIRREP)
       WRITE(6,1001)(VRTDOC(I),I=1,NIRREP)
       WRITE(6,1002)NAI
       WRITE(6,1003)N1I
       WRITE(6,1004)N2A
       WRITE(6,1005)NUMISCF
       WRITE(6,1006)NUMASCF
1001   FORMAT(8I5)
1002   FORMAT(T3,'Length of AI vector   : ',I5)
1003   FORMAT(T3,'Length of 1I vector   : ',I5)
1004   FORMAT(T3,'Length of 2A vector   : ',I5)
1005   FORMAT(T3,'Occupied RHF orbitals : ',I5)
1006   FORMAT(T3,'Virtual  RHF orbitals : ',I5)
      ENDIF
C 
C REORDER SCF EIGENVALUES BASED ON *RHF OCCUPATION*
C
      CALL GETREC(20,'JOBARC','QRHFEVAL',NBAS*IINTFP,SCR)
      IOFFIRR =1
      IOFFOCC2=NBAS+1
      IOFFVRT2=NBAS+NUMISCF+1
      DO 100 IRREP=1,NIRREP
       NOCC=POPRHF(IRREP)
       NVRT=VRTRHF(IRREP)
       IOFFOCC=IOFFIRR
       IOFFVRT=IOFFIRR+NOCC
       CALL SCOPY(NOCC,SCR(IOFFOCC),1,SCR(IOFFOCC2),1)
       CALL SCOPY(NVRT,SCR(IOFFVRT),1,SCR(IOFFVRT2),1)
       IOFFIRR=IOFFIRR+NOCC+NVRT
       IOFFOCC2=IOFFOCC2+NOCC
       IOFFVRT2=IOFFVRT2+NVRT
100   CONTINUE
      CALL PUTREC(20,'JOBARC','RHFEVAL ',NBAS*IINTFP,SCR(NBAS+1))
C      
      RETURN
      END
