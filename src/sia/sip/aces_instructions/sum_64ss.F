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
      subroutine sum_64ss(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c--------------------------------------------------------------------------
c   array1 --> 6D
c   array2 --> 4D
c   The last two indeces of array1 MUST be simple.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, k
      integer array, index, nindex, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk1, indblk2, get_block_index
      integer stack

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(6), findex(4)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer*8 addr, get_index_from_base
      double precision x(1)
      double precision y(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Detrmine the parameters of the first array: c_result_array
c---------------------------------------------------------------------------

      array = op(c_result_array)
      nindex = array_table(c_nindex, array)

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         sindex(i) = index
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, sop1(i), sop2(i))
      enddo

      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk1 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk1 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

c---------------------------------------------------------------------------
c   Detrmine the parameters of the second array: c_op1_array
c---------------------------------------------------------------------------

      array = op(c_op1_array)
      nindex = array_table(c_nindex, array)

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         findex(i) = index
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, fop1(i), fop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk2 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk2 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

c     write(6,*) ' ******************** '
c     write(6,*) ' OP1 :', (sindex(i), i=1,4)
c     write(6,*) ' OP2 :', (findex(i), i=1,4)

      call dosum_64ss(x(indblk1), sindex,
     *                  sop1(1), sop2(1), sop1(2), sop2(2),
     *                  sop1(3), sop2(3), sop1(4), sop2(4),
     *                  sop1(5), sop2(5), sop1(6), sop2(6),
     *                  x(indblk2), findex,
     *                  fop1(1), fop2(1), fop1(2), fop2(2),
     *                  fop1(3), fop2(3), fop1(4), fop2(4))

      return
      end

      subroutine dosum_64ss(x, sindex, a1, a2, b1, b2, c1, c2, d1, d2,
     *                                 s1, s2, ss1, ss2,
     *                        y, findex, e1, e2, f1, f2, g1, g2, h1, h2)
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer e1,e2,f1,f2,g1,g2,h1,h2
      integer s1, s2, ss1, ss2
      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(6), findex(4)

      integer i, j, k, a, b, c, d
      integer m, n, lda, ldb, ldc

      double precision x(a1:a2,b1:b2,c1:c2,d1:d2,s1:s2,ss1:ss2)
      double precision y(e1:e2,f1:f2,g1:g2,h1:h2), xtemp

      do d = d1, d2
      do c = c1, c2
      do b = b1, b2
      do a = a1, a2

c        xtemp = y(a,b,c,d)
c        do i = s1, s2
c        do j = ss1, ss2

            x(a,b,c,d,s1,ss1) = x(a,b,c,d,s1,ss1) +  y(a,b,c,d) ! xtemp
c           x(a,b,c,d,i,j) = x(a,b,c,d,i,j) +  xtemp

c        enddo
c        enddo

      enddo
      enddo
      enddo
      enddo

      return
      end

