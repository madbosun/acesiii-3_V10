      SUBROUTINE TOPT2S(ISPIN,ILIST,NTOP,NT2SIZ,NOCCA,NOCCB,
     &                  NVRTA,NVRTB,NSMSZ1,NSMSZ2,T,TOPT2,
     &                  SYVEC1,SYVEC2,ITOPT2,I,J,A,B,TYPE)
C
C THIS ROUTINE PICKS UP A T2 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT1 AND THEIR ASSOCIATED INDICES
C  IN THE I AND A VECTORS.
C
C INPUT:
C       ISPIN - THE LIST SUBTYPE
C       ILIST - THE LIST NUMBER
C       NTOP  - THE NUMBER OF LARGEST AMPLITUDES (BY ABSOLUTE VALUE)
C               WHICH ARE TO BE FOUND
C       NT2SIZ- THE TOTAL SIZE OF THE T2 VECTOR
C       NOCCA - THE NUMBER OF ALPHA OCCUPIED ORBITALS
C       NOCCB - THE NUMBER OF BETA OCCUPIED ORBITALS
C       NVRTA - THE NUMBER OF ALPHA VIRTUAL ORBITALS
C       NVRTB - THE NUMBER OF BETA VIRTUAL ORBITALS
C       NSMSZ1- THE TOTAL SIZE OF THE A,B SYMMETRY VECTOR
C               FOR THIS SPIN CASE (NVRTA*(NVRTA-1))/2 FOR
C               ISPIN=1, NVRTA*NVRTB FOR ISPIN=3, ETC.)
C       NSMSZ2- THE TOTAL SIZE OF THE I,J SYMMETRY VECTOR
C               FOR THIS SPIN CASE (NOCCA*(NOCCA-1))/2 FOR
C               ISPIN=1, NOCCA*NOCCB FOR ISPIN=3, ETC.)
C       SYVEC1- THE A,B SYMMETRY VECTOR
C       SYVEC2- THE I,J SYMMETRY VECTOR
C       TYPE  - ???????
C
C OUTPUT:
C       TOPT2 - THE NTOP LARGEST ELEMENTS IN THE T2 VECTOR
C               SORTED BY ABSOLUTE VALUE
C       I     - THE I INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       J     - THE J INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       A     - THE A INDICES CORRESPONDING TO THE VALUES IN TOPT1
C       B     - THE B INDICES CORRESPONDING TO THE VALUES IN TOPT1
C
C SCRATCH:
C
C       T     - USED TO HOLD THE SYMMETRY-PACKED T VECTOR
C       ITOPT2- USED TO HOLD THE OFFSETS CORRESPONDING TO THE
C               ELEMENTS IN TOPT2
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A(NTOP),B(NTOP),SYVEC1(NSMSZ1),SYVEC2(NSMSZ2)
cjp
      integer itrue, jtrue,atrue,btrue
      CHARACTER*2 SPCASE(3)
      CHARACTER*1 TYPE
      DIMENSION T(NT2SIZ),TOPT2(NTOP),ITOPT2(NTOP),I(NTOP),
     &          J(NTOP)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA SPCASE /'AA','BB','AB'/
cjp





cjp
cjp data for multireference state specific Brillouin-Wigner CC method
cjp coded by Jiri Pittner 1998-2000
cjp
      logical isbwcc,masik,isactive,bwgossip,useeq429,scfrefread
      logical bwwarning
      character*256 bwwarntext
      real*8 ecorrbw0,ecorrbw,epsilon0,fockcontr,denomblow,fockcd
      real*8 heff,heffevalr,heffevali,heffevecl,heffevecr,hdiagcontr
      real*8 fock2elcontr,enerscf,hcore,lambdahomotop,hfakt,diishonset
      real*8 fock2elcontr0,enerscf0,cbwstate,totmaxdenom,heffevecrold
      real*8 intruder,hfaktmax
      integer maxorb,maxref,nref, iref, iocc,iocc0
      integer iphnum,invpnum,invhnum,nbwstates,ibwstate
      integer ibwconvg,internfrom,internto,internnum,internindex
      integer ibwpass
      integer nactive,numactive,ihubaccorr,ihomotop
      integer internfrom1,internto1,internindex1,internnum1
      integer ihefferank, iheffefrom, iheffeto, iheffespin, maxexcit
      integer correctiontype
      integer maxbwwarnings, nproc, myproc
c
      parameter(maxorb=512,maxref=32,maxexcit=9,maxbwwarnings=10)
c      NOTE!!! change of maxorb parameter requires format change
c              and character* change in bwread routine!!!
      parameter(denomblow=1d250)
c

cjp common has been splitted in order to avoid problems
cjp with padding on different 32 and 64 bit architectures

      common/bwccint/isbwcc, masik, nref, iref,iocc(maxorb,maxref,2),
     +     iphnum(maxorb,maxref,2),invpnum(maxorb,maxref,2),
     +     invhnum(maxorb,maxref,2),
     +     isactive(maxorb,2),nbwstates,ibwstate(maxref+1),
     +     internfrom(maxref*(maxref-1)/2,maxref,3),
     +     internto(maxref*(maxref-1)/2,maxref,3),
     +     internindex(maxref*(maxref-1)/2,maxref,3),
     +     internnum(maxref,3),
     +     internfrom1(maxref*(maxref-1)/2,maxref,2),
     +     internto1(maxref*(maxref-1)/2,maxref,2),
     +     internindex1(maxref*(maxref-1)/2,maxref,2),
     +     internnum1(maxref,2),
     +     ibwpass,ibwconvg(maxref),bwgossip,useeq429,
     +     nactive(2),numactive(maxorb,2),ihubaccorr, ihomotop,
     +     iocc0(maxorb,2),scfrefread,
     +     ihefferank(maxref,maxref),iheffefrom(maxexcit,maxref,maxref),
     +     iheffeto(maxexcit,maxref,maxref),
     +     iheffespin(maxexcit,maxref,maxref),
     +     correctiontype,bwwarning(maxbwwarnings),
     +     bwwarntext(maxbwwarnings),nproc,myproc

      common/bwccreal/ecorrbw,epsilon0,cbwstate(maxref+1),
     +     fockcontr(maxorb*(maxorb+1)/2,2),fockcd(maxorb,maxref,2),
     +     heff(maxref,maxref),heffevalr(maxref),heffevali(maxref),
     +     heffevecl(maxref,maxref),heffevecr(maxref,maxref),
     +     hdiagcontr(maxref), fock2elcontr(maxorb,2),
     +     enerscf(maxref), hcore(maxorb,2),
     +     lambdahomotop,hfakt,diishonset,enerscf0,
     +     fock2elcontr0(maxorb,2),ecorrbw0,totmaxdenom,
     +     heffevecrold(maxref,maxref),intruder,hfaktmax

c
c
cjp BRIEF DESCRIPTION OF VARIABLES INTRODUCED FOR THE MR-BWCC ROUTINES
cjp IN FACT, A LOT OF THAT COULD BE USEFUL FOR ANY HILBERT-SPACE MR-CC
c
c
c isbwcc ... flag for doing bwcc calculation
c maxbwwarnings, bwwarning, bwwarntext ... serious warnings will be
c    summarized at the end of xvcc output for the user's' convenience
c ihefferank(jref,iref) ... degree of excitation between jref and iref
c iheffefrom(maxexcit,jref,iref) , iheffeto, iheffespin ... list of
c    indices of that excitation, sorted according to spin and then the
c    indices, numbers stored are defined as effective particle-hole
c    indices of reference iref
c ihomotop ... whether to use homotopic transition to the
c    size-extensivity correction, after which iteration (if .ne.0)
c lambdahomotop scaling factor of the geometrical series of
c    lambda 1->0 transition
c hfakt ... current value of the homotopy parameter
c hfaktmax ... maximal homotopy parameter allowed to consider cc
c    equations converged
c diishonset ... at which value of hfact restart diis convergence acceleration
c masik ... prepare sorted integral file for the program by Masik and stop
c nref ... number of reference configurations
c iref ... current reference configuration and fermi vacuum
c bwgossip ... switch on debugging output
c ibwpass ... routines like newt2 have to be splitted in two passes -
c    construction of Heff and amplitude update after heff is diagonalized
c    for backw. compatibility, instead of introducing a new routine
c    the same routine does different things being called twice with
c    different ibwpass value
c ecorrbw ... correlation energy from BWCC - Heff(iref,iref) ...
c    denominator correction
c ecorrbw0 .... ecorrbw, but not scaled by the homotopic factor hfact
c denomblow ... huge number to cause division underflow - used for
c    zeroing out the internal amplitudes automatically
c nactive(spin): total count of active spinorbitals
c numactive(i=1..nactive,spin): number of i-th active spinorbital
c    in sequential numbering
c isactive(maxorb,spin): belongs given orbital to the active space?
c    for RHF, the beta ones must be initialized to be identical with alpha ones
c iocc(maxorb,1..nref,spin): defines the nref reference configurations
c    for both RHF and UHF:  iocc(i,iref,spin)=0 or 1
c iphnum(orbital no, iref, spin): gives the effective number of orbital
c    (both particle and hole ones are counted starting from 1)
c invpnum(eff.p.orb.,iref,spin): gives true orbital no. from the
c    effective particle one
c invhnum(eff.h.orb.,iref,spin): gives true orbital no. from the
c    effective hole one
c    all these three ones must be in RHF case initialised to be equal in the
c    alpha and beta parts to keep the code unique and simple
c internfrom(sequence counter n,iref,ispin) is for the ab spin case, the
c    other ones have to be iuhf-indexed
c    internfrom, internto: they are first and second index of n-th internal
c    excitation when processing reference given by second index to the array
c
c internindex(sequence counter,iref,ispin) is the position of
c    corresponding denominator in the denominator list
c internnum(iref,ispin) - number of internal excitation in that
c    category = max sequence counter here ispin=1,2,3 for AA,BB,AB
c internfrom1 etc. are analogous quantities for monoexcitations, here
c    ispin =1,2 note for later: all intern.... quantities are irrep-specific!
c fockcontr(findex(i,j),ispin) ... addition to the fock matrix of
c    reference no.1 to obtain the fock matrix of current reference
c    (fermi vacuum)
c fockcd(i,iref,ispin) ... diagonal part of that correction for ref. no. iref.
c hcore(i,ispin) ... one electron diagonal hamiltonian elements
c fock2elcontr(i,ispin) ... 2el contribution to the diagonal fock element
c    used temporarily
c hdiagcontr(iref) ... contribution of differences of HF energies of
c    different Fermi vacua to diagonal Heff elements
c enerscf(iref) ... HF energy of iref-th fermi vacuum
c ihubaccorr ... =1 ... calculate the size extenzivity correction for BWCC
c                =2,3 ... second and third pass of that calculation
c iocc0(maxorb,2) ... like iocc, but for dummy reference configuration
c    corresponding to SCF WF
c fock2elcontr0(maxorb,2) ... like fock2elcontr0 but for dummy reference
c    configuration
c enerscf0 ... like enerscf, but for dummy reference config
c scfrefread ... tells to bwprep that SCF reference has been read from input
c    and should not be generated automatically from nocc(ispin)
c nbwstates ... how many states to average
c ibwstates(1..nbwstates),cbwstates() ... their numbers and coefficients
c correctiontype ... 0=DC,L T2 term is removed/scaled, 1=DC/L term is not
c    removed/scaled
c totmaxdenom ... max 1/denom found for given reference's' fermi vacuum -
c    as indication of possible intruder problem
c intruder ... limit of 1/denom to be considered intruder and its
c    amplitude zeroed
c for parallelization
c nproc ... number of processors (counted from 1)
c myproc ... number of the processor currently executing the code
c    (counted from 1)



C
C THIS ROUTINE PICKS UP A T2 VECTOR FROM MOIO(ISPIN,ILIST) AND RETURNS
C  THE NTOP LARGEST ELEMENTS IN TOPT2 AND THEIR ASSOCIATED INDICES
C  IN THE I,J,A AND B VECTORS.
C
C
C STATEMENT FUNCTIONS
C
      ILRG(INTX)=INT(0.5D0*(1.D0+SQRT(8.0*INTX-7)))+1
      NNM1O2(IX)=(IX*(IX-1))/2
C
      IF(NT2SIZ.EQ.0)RETURN
C
      CALL GETALL(T,NT2SIZ,1,ILIST)
C
C COMPUTE NORM OF T2 VECTOR
C
      X=SNRM2(NT2SIZ,T,1)
      T2NORM=X
C
C PUT NTOP LARGEST ELEMENTS INTO TOPT2 AND THEIR ASSOCIATED INDICES INTO
C   ITOPT2.
C
      CALL SCANVC(T,TOPT2,ITOPT2,NTOP,NT2SIZ)
C
C BRANCH HERE FOR DIFFERENT SPIN CASES, BECAUSE INDICES HAVE TO BE
C   FIRURED OUT FROM ITOPT2.
C
      DO 10 IRANK=1,NTOP
       CALL PCKIND(ITOPT2(IRANK),ISPIN,SYVEC1,SYVEC2,NSMSZ1,NSMSZ2,
     &             I(IRANK),J(IRANK),A(IRANK),B(IRANK))
10    CONTINUE
cjp
      if(isbwcc) then
cjp prepare spin indices according to ispin case
        if(ispin.eq.1) then
          isp1=1
          isp2=1
        else
          if(ispin.eq.2) then
            isp1=2
            isp2=2
          else
            isp1=1
            isp2=2
          endif
        endif
        write(*,202)type,SPCASE(ISPIN),iref
        else
      WRITE(*,200)TYPE,SPCASE(ISPIN)
        endif
200   FORMAT(T3,' Largest ',A1,'2 amplitudes for spin case ',A2,':')
202   format(T3,' Largest ',A1,'2 amplitudes for spin case ',A2,
     &    ' of reference no.: ',i3)
      IF(ISPIN.EQ.3)WRITE(*,299)
      WRITE(*,300)
299   FORMAT(3X,3(' ',3X,'_',3X,' ',3X,'_',13X))
300   FORMAT(3X,3('i',3X,'j',3X,'a',3X,'b',13X))
      WRITE(*,400)
400   FORMAT(77('-'))
cjp
      if(isbwcc) then
cjp ... this is not OK for uhf-cc, ignore for now
       if(bwgossip) then
      write(*,502) (irank,i(irank),j(irank),a(irank),b(irank),
     &   invhnum(i(irank),iref,isp1),
     & invhnum(j(irank),iref,isp2),invpnum(a(irank)-nocca,iref,isp1),
     & invpnum(b(irank)-nocca,iref,isp2),topt2(irank),irank=1,ntop)
        else
        write(*,500) (invhnum(i(irank),iref,isp1),invhnum(j(irank),
     &  iref,isp2),invpnum(a(irank)-nocca,iref,isp1),
     & invpnum(b(irank)-nocca,iref,isp2),topt2(irank),irank=1,ntop)
        endif
      else
      WRITE(*,500)(I(IRANK),J(IRANK),A(IRANK),B(IRANK),
     &                 TOPT2(IRANK),IRANK=1,NTOP)
      endif
      WRITE(*,400)
      WRITE(*,600)TYPE,SPCASE(ISPIN),NT2SIZ,T2NORM
      WRITE(*,400)
cpn
      if(isbwcc.and.bwgossip) then
c This will print all the amplitudes and their indexes in order
      write(*,203)type,spcase(ispin),iref
 203  format(t3,' List of all ',A1,'2 amplitudes for spin case ',A2,
     &    ' of reference no.: ',i3)
      if(ispin.eq.3)write(*,299)
      write(*,300)
      write(*,400)
      do 11 irank = 1, nt2siz
         call pckind(irank,ispin,syvec1,syvec2,nsmsz1,nsmsz2,
     &        i(1),j(1),a(1),b(1))
         itrue=invhnum(i(1),iref,isp1)
         jtrue=invhnum(j(1),iref,isp2)
cjp ... this is not OK for uhf-cc, ignore for now
         atrue=invpnum(a(1)-nocca,iref,isp1)
         btrue=invpnum(b(1)-nocca,iref,isp2)
         write(*,502)irank,i(1),j(1),a(1),b(1),itrue,jtrue,atrue,
     &         btrue,t(irank)
 11   continue
 502  format(3x,i5,6x,4i5,': true ind : ',4i5,':',f12.8)
      endif
cpn end
600   FORMAT(T3,' Norm of ',A1,'2',A2,' vector (',I9,
     &       ' symmetry allowed elements):',F14.10,'.')
500   FORMAT((2('[',I3,1X,I3,1X,I3,1X,I3,']',F8.5,1X),'[',I3,1X,
     &I3,1X,I3,1X,I3,']',F8.5))
      RETURN
      END
