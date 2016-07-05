C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine clear_saved_requests()
c-------------------------------------------------------------------------
c   Clear the saved requests that were stored to this point.
c-------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'saved_data.h'
      include 'int_gen_parms.h'
      integer status(MPI_STATUS_SIZE)
      integer i, ierr

      common /saved_requests/my_saved_requests(max_procsx)
      integer my_saved_requests

      do i = 1, n_saved_requests
         call mpi_wait(my_saved_requests(i), status, ierr)
      enddo

      n_saved_requests = 0
      return
      end
