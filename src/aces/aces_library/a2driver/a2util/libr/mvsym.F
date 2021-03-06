
C MOVES SYMMETRY PACKED LIST LIST1 TO SYMMETRY PACKED LIST LIST2.

      SUBROUTINE MVSYM(ICORE,MAXCOR,IUHF,ICACHE,LIST1,LIST2)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),IJUNK(18)
      NSIZE = 0
      ILIST1 = LIST1
cYAU      IF (LIST1.GT.100) ILIST1 = ILIST1
      DO IRREP = 1, NIRREP
         NUMSYM = IRPDPD(IRREP,ISYTYP(2,ILIST1))
         DSZSYM = IRPDPD(IRREP,ISYTYP(1,ILIST1))
         NSIZE  = NSIZE + ( DSZSYM * NUMSYM )
      END DO
      CALL GETALL(ICORE,NSIZE,ICACHE,LIST1)
      CALL PUTALL(ICORE,NSIZE,ICACHE,LIST2)
      RETURN
      END
