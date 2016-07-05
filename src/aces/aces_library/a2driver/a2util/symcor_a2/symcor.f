
C THIS PROGRAM DETERMINES SYMMETRY ADAPTED COORDINATES WHICH
C ARE USED FOR NUMERICAL HARMONIC FREQUENCY CALCULATIONS.  IT
C CAN EASILY BE ADAPTED TO CONSTRUCT OTHER SYMMETRY ADAPTED
C QUANTITIES, SUCH AS BASIS FUNCTIONS.  IT CURRENTLY WORKS FOR
C ALL POINT GROUPS OTHER THAN THOSE WHICH CONTAIN COMPLEX
C REPRESENTATIONS (CN, SN, CNH (N>2) AND T AND TH).  FOR THESE
C "DANGEROUS GROUPS", THE ABELIAN SUBGROUP IS USED.
C WRITTEN BY J.F. STANTON, 1991-1992.

c The entire operating procedure of symcor is explained in README.




































































































































































































































































































































































































































































































































      SUBROUTINE SYMCOR(ICORE,ICRSIZ)
      IMPLICIT INTEGER (A-Z)
      INTEGER ICORE(ICRSIZ), I0

C The dimension of IRRNM array is changed to 14 to accomodate the
C largest possible dimension. Note that in joda D_infh and C_infv
C are set to D_8h and C_8v. Ajith 08/2001

      DOUBLE PRECISION E, TRAD, PRD_ENRG_CHNG, TAU, TotEng
      CHARACTER*8 GRPSYM,IRRNM(14)
      CHARACTER*4 DOIT
      LOGICAL FIRST_PASS,GRAD_EXIST,I_HAVEAGEOM

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100), IFLAGS2(500)
c      COMMON /FLAGS2/ IFLAGS2(500)
      LOGICAL          ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT,
     &                 SNPTGRAD
      COMMON /CONTROL/ ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT,
     &                 SNPTGRAD 
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end

      DATA IZERO /0/

c   o initialize the free-space index in iCore
      I0=1

C The values get assigned for PASS1 JOBARC record has very important
C function in joda when it is executing any calculation that involve
C numerical derivatives.
C
C When the execution transfers to here for the first time, joda
C has run one time as a single point calculation. There are no
C JOBARC records GRADIENT or PASS1. So GRAD_EXIST is false,
C PASS1 is 0 and GMTRYOPT is true or false depending on whether
C we are doing geometry optimizations or vibrational frequencies.
C
C Numerical geometry optimization - PASS1 is set to 1 until
C full set of gradients available. Then it is set to 0.
C
C Analytical gradient frequencies  - PASS1 is set to 1 in the
C very first time. It remains as 1 until the full Hessian is
C is available, then it is set to 0.
C
C Numerical gradient frequencies  - PASS1 get set to -1 in the
C first run. Before the next execution joda, this  get reset
C back to 1. When the full Hessian is available it get reset to 0.
C (-1 in the first call tell the code that it has to do energy
C points). 03/2003 Ajith Perera.

      CALL GETREC(-1,'JOBARC','PASS1',1,PASS1)
      CALL GETREC(0,'JOBARC','GRADIENT',LENGTH,ITMP)
      Call GETREC(-1,  'JOBARC', 'HAVEGEOM',1 , LEN_HAVEGEOM)
      I_HAVEAGEOM = (LEN_HAVEGEOM .EQ. 1)
      GRAD_EXIST  = LENGTH.NE.-1

c   o initialize /CONTROL/
      ENERONLY = IFLAGS2(138).EQ.2
      GRADONLY = IFLAGS2(138).EQ.1
      ROTPROJ  = IFLAGS (80).EQ.0
      RAMAN    = IFLAGS2(151).EQ.1
      GMTRYOPT = IFLAGS2(5).NE.0 .AND. .NOT.
     &           I_HAVEAGEOM
      SNPTGRAD = .NOT. (GMTRYOPT.OR.IFLAGS(54).NE.0)

c   o detect dangerous point groups and avoid them like hell
      CALL GETCREC(20,'JOBARC','FULLPTGP',8,GRPSYM)
      IPOS = LINBLNK(GRPSYM(1:4))
      IF (PASS1.NE.1.AND.
     &    (GRPSYM(1:1).EQ.'C'.AND.GRPSYM(IPOS:IPOS).NE.'v'.OR.
     &     GRPSYM(1:3).EQ.'T h'.OR.
     &     GRPSYM(1:3).EQ.'T  '.OR.
     &     GRPSYM(1:1).EQ.'S'
     &    )
     &   ) THEN
         DOIT='COMP'
         WRITE(6,410)
 410     FORMAT(T3,'@SYMCOR: Full point group is dangerous! ',
     &             'Abelian subgroup will be used.')
         ITMP=1
         CALL PUTREC(20,'JOBARC','DANGERUS',1,ITMP)
      ELSE
         CALL GETREC(-1,'JOBARC','DANGERUS',1,ITMP)
         IF (IFLAGS(79).EQ.0.AND.ITMP.EQ.0) THEN
            DOIT='FULL'
         ELSE
            DOIT='COMP'
            ITMP=0
         END IF
         IF (PASS1.NE.1) CALL PUTREC(20,'JOBARC','DANGERUS',1,ITMP)
      END IF

      CALL GETREC(20,'JOBARC','NATOMS',1,NATOM)

      IF (PASS1.EQ.-1) THEN
         IF (IFLAGS(87).EQ.0) THEN
            CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,E)
         ELSE
            CALL GETREC(20,'JOBARC','TOTENER2',IINTFP,E)
         END IF
         CALL PUTREC(20,'JOBARC','REFENERG',IINTFP,E)
         PASS1=1
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
      ELSE
         E = 0.d0
      END IF

      FIRST_PASS = ( (PASS1.EQ.0) .OR. (GRAD_EXIST.AND.GMTRYOPT) )
      Print*, "Geo opt flag set?, vlaue of pass1, gradient exist",
     &        "geometry optimization?, energy only?:"
      Print*, IFLAGS2(5), PASS1, GRAD_EXIST, GMTRYOPT, 
     &        ENERONLY
      IF (FIRST_PASS) THEN

c      o zero icore so we can use it to initialize some records
         do i = 1, iintfp*3*natom
            icore(i) = 0
         end do

c      o create a record for the force constants so we can use them in testing
         IF (.NOT.GMTRYOPT) THEN
            CALL PUTREC(20,'JOBARC','FORCECON',IINTFP*3*NATOM,ICORE)
         END IF

c      o create records that must survive truncating JOBARC
         CALL PUTREC(20,'JOBARC','REFENERG',IINTFP,E)
         CALL PUTREC(20,'JOBARC','LASTGEOM',1,IZERO)
c      - these records are used by the new geometry stepping algorithms
         CALL GETREC(-1,'JOBARC','T_RADIUS',IINTFP,TRAD)
         CALL GETREC(-1,'JOBARC','PRDENCHN',IINTFP,PRD_ENRG_CHNG)
         CALL GETREC(-1,'JOBARC','PRVIUTAU',IINTFP,TAU)
         CALL GETREC(-1,'JOBARC','OLDENERG',IINTFP,TotEng)
         CALL PUTREC(20,'JOBARC','T_RADIUS',IINTFP,TRAD)
         CALL PUTREC(20,'JOBARC','PRDENCHN',IINTFP,PRD_ENRG_CHNG)
         CALL PUTREC(20,'JOBARC','PRVIUTAU',IINTFP,TAU)
         CALL PUTREC(20,'JOBARC','OLDENERG',IINTFP,TotEng)
c         - these three records are used by vee and vcceh
            CALL PUTREC(20,'JOBARC','PRINSPIN',1,IZERO)
            CALL PUTREC(20,'JOBARC','PRINFROM',IINTFP,E)
            CALL PUTREC(20,'JOBARC','PRININTO',IINTFP,E)
c      - these are used for "restarts"
         CALL PUTCREC(20,'JOBARC','PREVPTGP',4,'NONE')
         CALL PUTREC(0,'JOBARC','TGSSOCCA',8,0)
         CALL PUTREC(0,'JOBARC','TGSSOCCB',8,0)

c      o initialize the finite difference grid
         call init_fd(doit,natom,icore,icrsiz,nener,
     &                IFLAGS(1).GE.10)

         IF (.NOT.GMTRYOPT.AND.NENER.NE.0.AND.IRANK.EQ.0) THEN
c         o do the reference geom iff this is a vib freq calc with numerical
c           gradients and it is the root process
            PASS1=-1
         ELSE
            PASS1=1
         END IF
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
C The following record no matter how awkward it looks is essential. Note at
C this point we ran xjoda and there are certain records in JOBARC and then we
C add the "NEXTGEOM" record as a marker (initialized to junk at icore). Then we
C continue with vmol2ja and so on. In the next pass we are going  to fill the
C "NEXTGEOM" (see below the call to NEXTGEO) with proper data. The call to
C ACES_JA_TRUNCATE would delete all the records that comes after "NEXTGEOM"
C record. If we don't fixed the position of the "NEXTGEOM" record after the
C first JODA run, the records that are dependent on the symmetry get  written
C to the JOBARC before the "NEXTGEOM" appear and ACES_JA_TRUNCATE call would do
C nothing to erase them. Then When it times to run the lower symmtry points,
C the PUTREC is going to cry since we will be trying to change the record
C lengths!  Ajith Perera 08/2000.
         CALL PUTREC(20,'JOBARC','NEXTGEOM',3*NATOM*IINTFP,ICORE)
         IF (PASS1.EQ.-1) RETURN

      END IF

C THIS MUST BE CALLED ON EVERY PASS
C-----------------------------------------------------------------------
CJDW May/June 1996. Modified allocation for subroutine UPD_FD.
CJDW/AP Jul/1998. Extensions for Raman intensities.
C-----------------------------------------------------------------------

      I000=I0
      I010=I000 + IINTFP*9*NATOM*NATOM
      I020=I010 + IINTFP*9*NATOM*NATOM*3*NATOM
      I030=I020 + IINTFP*9*NATOM*NATOM*3
      I040=I030 + IINTFP*9*NATOM*NATOM*9
      I200=I040 + 9*NATOM*NATOM + IAND(NATOM,1)

      CALL UPD_FD(NATOM,DOIT,IMORE,
     &            ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &            ICORE(I040),ICORE(I200),(ICRSIZ-I200+I0)/IINTFP)

      IF (IMORE.EQ.1) THEN
c      o if there are more points, delete all the geometry-dependent records
         CALL ACES_JA_TRUNCATE('NEXTGEOM',1)
      ELSE
         CALL GETREC(-1,'JOBARC',DOIT//'NIRX',1,NIRREPF)
         IF (GMTRYOPT .OR. SNPTGRAD) THEN
            I010=I0
            I020=I010 + 3*NATOM + IAND(NATOM,1)
            I030=I020 + IINTFP*3*NATOM
            I200=I030 + IINTFP*3*NATOM
            CALL SETGRD(NATOM,NIRREPF,DOIT,
     &                  IRRNM,ICORE(I010),
     &                  ICORE(I020),ICORE(I030),
     &                  ICORE(I200),(ICRSIZ-I200+I0)/IINTFP
     &                 )
         ELSE
            CALL GETREC(20,'JOBARC','ORDERREF',1,IORDERF)
            I010=I0
            I011=I010 + 3*NATOM
            I012=I011 + 3*NATOM
            I020=I012 + IINTFP*9*IORDERF
            I030=I020 + IINTFP*9*NATOM*NATOM
            I040=I030 + IINTFP*9*NATOM*NATOM
            I050=I040 + IINTFP*9*NATOM
            I060=I050 + IINTFP*9*NATOM
            I070=I060 + IINTFP*9*NATOM*3
            I200=I070 + IINTFP*9*NATOM*3
            CALL SETFCM(NATOM,NIRREPF,IORDERF,DOIT,
     &                  IRRNM,ICORE(I010),ICORE(I011),ICORE(I012),
     &                  ICORE(I020),ICORE(I030),
     &                  ICORE(I040),ICORE(I050),
     &                  ICORE(I060),ICORE(I070),
     &                  ICORE(I200),(ICRSIZ-I200+I0)/IINTFP
     &                 )
         END IF

c      o signal the end of this findiff series
         call putrec(1,'JOBARC','FNDFDONE',1,1)

         PASS1=0
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
      END IF

c   o mark the symmetry of the previous geometry
      CALL GETCREC(20,'JOBARC','COMPPTGP',4,DOIT)
      CALL PUTCREC(20,'JOBARC','PREVPTGP',4,DOIT)

      RETURN
      END

