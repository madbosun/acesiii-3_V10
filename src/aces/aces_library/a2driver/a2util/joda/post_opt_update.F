      Subroutine Post_opt_update
c
      Implicit Double Precision (A-H, O-Z)
c
#include "mxatms.par"
#include "flags.h"
#include "baslen.par"
#include "jodaflags.com"
#include "machsp.com"
#include "coord.com"
c
      Integer Genby(Mxatms) 
      Double Precision Orient(3, 3)
      Character*(baslen) BasNam(Mxatms)
      Character*4 FPGrp, BPGrp, PGrp 
      Character*5, Zsym(Mxatms)
      Dimension Scratch(Mxatms) 
      Logical Can_do_freq 
c
      Common /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
c
      Call Entry(BasNam, .True., Can_do_freq)
c
      Call Getrec(20, 'JOBARC', 'ZMATATMS', 1, NATOMS)
      Call Getrec(20, 'JOBARC', 'LINEAR  ', 1, ILINEAR)
      If (ILINEAR .EQ. 1) Then
          NX   = 3*NATOMS
          NXM6 = NX - 5
      Else
          NX   = 3*NATOMS
          NXM6 = NX -6
      Endif  
c
      Call Getrec(20, "JOBARC", "COORD   ", NX*IINTFP, Q)
      Call Getrec(20, 'JOBARC', 'CONCTVTY', NX, NCON)
      Call Getrec(-1, 'JOBARC', 'CORD_INT', NX*IINTFP, R)
      Call Bohr2angs(R, NX)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', NATOMS*IINTFP,
     &           ATMASS)
      Call Getrec(20, 'JOBARC', 'ORIENTMT', 9*IINTFP,
     &            ORIENT)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', NATOMS, IATNUM)
      Call Getrec(20, 'JOBARC', 'ICSQUASH', NX, ISQUASH)
      Call Getcrec(20, 'JOBARC', "PTGP    ", 4, FPGRP)
      Call Getcrec(20, 'JOBARC', "ABL_PTGP", 4, BPGRP)
      Call Getcrec(20, 'JOBARC', "CMP_PTGP", 4, PGRP)
      Call Getcrec(20, 'JOBARC', "ZSYM", 5*NATOMS, ZSYM)
#ifdef _DEBUG_LVL0
           Print*, "The data read from JOBARC in post_opt_update"
           Print*, "The NATOMS:", NATOMS
           Print*, "Internal coords:",(R(I), I=1, NX)
           Print*, "The Cartesian coords:", (Q(I), I=1, NX) 
           print*, "The connectivities:", (NCON(I), I=1, NX)
           Print*, "The atomic charges:", (IATNUM(I), I=1, NX/3)
           Print*, "The point group", PGRP 
#endif
c
c The scratch and Nunique are not used. Genby, istat and Genby
c are internal 

      Call MkVMOL(Q, PGrp, NAtoms, NUnique, ZSym, IAtNum,
     &            GenBy, Scratch, IStat, BasNam)
      
#ifdef _DEBUG_LVL0
      Print*, "out from Mkvmol in post_opt_update"
#endif
      Return
      End


