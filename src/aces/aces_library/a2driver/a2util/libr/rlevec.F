      SUBROUTINE RLEVEC(ICORE,MAXCOR,IUHF,TAU,IORDER,SINGLE,E,
     &                  RLETYP,MAXORD)
C
C THIS ROUTINE FORMS AN EXTRAPOLATED T VECTOR ACCORDING TO THE
C  RLE TAU VECTOR.
C
C   T2(new) = [1 - SUM TAU(k)] T(1) + SUM TAU(k) * T(k)
C                   k                  k
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION Z,COEFF,TAU,E,SDOT,Z0
      LOGICAL SINGLE,NONHF
      DIMENSION ICORE(MAXCOR),TAU(IORDER)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /RLEIO/ IOLDEST
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      NONHF=.FALSE.
      IF(IFLAGS(32).NE.0.OR.IFLAGS(33).NE.0.OR.IFLAGS(38).NE.0)THEN
       NONHF=.TRUE.
      ENDIF
C
C FORM  Z= 1 - SUM TAU(k)
C                k
C
      Z=1.0
      DO 10 I=1,IORDER
       Z=Z-TAU(I)
10    CONTINUE
      IF(RLETYP.EQ.2)E=0.0
      DO 100 ISPIN=3,1,IUHF-2
       SIZET=ISYMSZ(ISYTYP(1,60+ISPIN),ISYTYP(2,60+ISPIN))
       IF(SINGLE.AND.ISPIN.LT.3)SIZET=SIZET+NT(ISPIN)
       I000=1
       I010=I000+IINTFP*SIZET
       I020=I010+IINTFP*SIZET
       NLIST=60+ISPIN
       CALL IZERO(ICORE(I000),SIZET*IINTFP)
       IF(MOD(RLETYP,2).EQ.0)THEN
C
C FORM  Z*T(1)
C 
        CALL GETLST(ICORE(I000),1,1,1,ISPIN,89)
        CALL SSCAL(SIZET,Z,ICORE(I000),1)
C
C NOW PICK UP ALL OTHER T VECTORS, WEIGHT THEM WITH TAU AND
C  ACCUMULATE INTO BOTTOM OF CORE.
C
       ENDIF
       DO 20 I=1,IORDER
        COEFF=TAU(I)
        CALL GETLST(ICORE(I010),IGET(I),1,1,ISPIN,70)
        CALL SAXPY(SIZET,COEFF,ICORE(I010),1,ICORE(I000),1)
20     CONTINUE
       IF(RLETYP.NE.2.AND.RLETYP.NE.3)THEN
        CALL PUTLST(ICORE(I000),IGET(1),1,1,ISPIN,70)
       ELSE
        CALL PUTLST(ICORE(I000),IGET(1),1,1,ISPIN+3,70)
       ENDIF
C
C NOW OVERWRITE LISTS 61-63 WITH THE UPDATED T2 VECTOR, AND THE
C   T1 LISTS IF THE MODEL IS NOT CCD OR LCCD.
C
       IOFF=1
       DO 30 IRREP=1,NIRREP
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,NLIST))
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,NLIST))
        CALL PUTLST(ICORE(IOFF),1,NUMDIS,1,IRREP,NLIST)
        IOFF=IOFF+NUMDIS*DISSIZ*IINTFP
30     CONTINUE
       IF(SINGLE.AND.ISPIN.LT.3)THEN
        CALL PUTLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
       ENDIF
C
C COMPUTE T1+T2 ENERGY FOR RLE/NOJACOBI METHOD
C
       IF(RLETYP.EQ.2)THEN
        LENT2=(IOFF-1)/IINTFP
        LENT1=SIZET-LENT2
        CALL GETALL(ICORE(I010),LENT2,1,13+ISPIN)
        Z0=SDOT(LENT2,ICORE(I000),1,ICORE(I010),1)
        IF(IUHF.EQ.0.AND.ISPIN.EQ.1)Z0=2.0*Z0
        E=E+Z0
        IF(NONHF.AND.SINGLE.AND.ISPIN.LT.3)THEN
         CALL GETLST(ICORE(I010),1,1,1,2+ISPIN,93)
         Z0=SDOT(LENT1,ICORE(IOFF),1,ICORE(I010),1)
         IF(IUHF.EQ.0.AND.ISPIN.EQ.1)Z0=2.0*Z0
         E=E+Z0
        ENDIF
       ENDIF
100   CONTINUE
      RETURN
      END 
