
c This routine reads a physical record from a direct access file.

c INPUT
c int IUNIT   : the external unit number of the direct access file
c int IREC    : the index of the physical record to read
c int ILENGTH : the integer-length of the record and destination array

c OUTPUT
c int IDEST : the contents of the physical record

      subroutine aces_io_read(iUnit,iRec,iDest,iLength)
      implicit none

c ARGUMENTS
      integer iUnit, iRec, iLength
      integer iDest(iLength)

c INTERNAL VARIABLES
      integer iStat

c ----------------------------------------------------------------------


      iStat = 0
c   o assert iUnit and iRec are natural
      if ((iUnit.lt.1).or.(iRec.lt.1)) then
         print *, '@ACES_IO_READ: Assertion failed.'
         print *, '   iUnit = ',iUnit
         print *, '   iRec  = ',iRec
         iStat = 1
      end if
c   o assert iLength is whole
      if (iLength.lt.0) then
         print *, '@ACES_IO_READ: Assertion failed.'
         print *, '   iLength = ',iLength
         iStat = 1
      end if
      if (iStat.ne.0) call aces_exit(iStat)


c ----------------------------------------------------------------------

      read(unit=iUnit,rec=iRec,err=666,iostat=iStat) iDest
      return

 666  call aces_io_error('ACES_IO_READ',iUnit,iStat)

c     end subroutine aces_io_read
      end

