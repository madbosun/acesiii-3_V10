
C THIS ROUTINE WRITES OUT INFORMATION WHICH ALLOWS LATER MODULES TO
C RUN SUCCESSFULLY. THIS MUST BE CALLED BY THE LAST EXECUTABLE STATEMENT
C IN ALL PROGRAMS.

      SUBROUTINE CRAPSO
      IMPLICIT NONE

      CHARACTER*80 FNAME
      INTEGER ILENGTH
      LOGICAL bExist,bOpened

      CALL GFNAME('MOINTS  ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=bExist,OPENED=bOpened)
      IF (bExist.AND.bOpened) CALL TERMIO
      CALL ACES_JA_FIN

      RETURN
      END

