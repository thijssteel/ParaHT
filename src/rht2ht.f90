module rht2ht
   use omp_lib
   use utils

contains

   subroutine rht2ht_householder_nowy( n, A, lda, B, ldb, Q, ldq, Z, ldz, r, nq, time_components )
      implicit none
      ! arguments
      integer, intent(in) :: n, lda, ldb, ldq, ldz, r, nq
      double precision, intent(inout) :: A(lda,*), B(ldb,*), Z(ldz,*),Q(ldq,*)
      logical :: time_components

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, parameter :: min_blocksize = 32

      ! local variables
      integer :: i, jc, jb, ierr, nb, nblocks, jblock, ic, ib, istart, istop, istart2,&
         istop2, r1, r2, c1, c2, nthreads, blocksize
      double precision, allocatable :: U(:,:,:), V(:,:,:), &
         work(:), B_copy(:,:), tauU(:,:), tauV(:,:),&
         U2(:,:,:), V2(:,:,:),tauU2(:,:), tauV2(:,:)

      ! temporary timing variables
      logical, parameter :: time_tasks = .false.
      integer(8) :: starttime,endtime,count_rate
      double precision :: generationtime, lookaheadtimeA, lookaheadtimeB,&
         updatetimeQ, updatetimeA, updatetimeB, updatetimeZ
      double precision :: generationtime2, applicationtime, lookaheadtime

      ! ---------------
      ! executable code
      ! ---------------

      nblocks = (n-1)/r
      if( mod( (n-1),r ) .ne. 0 ) nblocks = nblocks + 1


      allocate(work( (r+nq)*n ))
      allocate(B_copy( r, r ))
      allocate( U( r+nq, nq, nblocks ) )
      allocate( V( r+nq, nq, nblocks ) )
      allocate( tauU( nq, nblocks ) )
      allocate( tauV( nq, nblocks ) )
      allocate( U2( r+nq, nq, nblocks ) )
      allocate( V2( r+nq, nq, nblocks ) )
      allocate( tauU2( nq, nblocks ) )
      allocate( tauV2( nq, nblocks ) )

      !$omp parallel
      !$omp single

      generationtime = 0.0d0
      generationtime2 = 0.0d0
      lookaheadtime = 0.0d0
      applicationtime = 0.0d0
      lookaheadtimeA = 0.0d0
      lookaheadtimeB = 0.0d0
      updatetimeA = 0.0d0
      updatetimeB = 0.0d0
      updatetimeQ = 0.0d0
      updatetimeZ = 0.0d0

      nthreads = omp_get_max_threads()
      blocksize = n/nthreads
      if( mod( blocksize, 16 ) .ne. 0 ) blocksize = blocksize + (16-mod( blocksize, 16 ))
      blocksize = max( min_blocksize, blocksize )
      blocksize = min( n, blocksize )

      ! jblock indicates the first column of the block we are reducing
      do jblock = 1,n-2,nq
         ! do jblock = 1,1,nq

         nblocks = (n-jblock)/r
         if( mod( (n-jblock),r ) .ne. 0 ) nblocks = nblocks + 1

         if( time_components ) then
            !$omp taskwait
            call system_clock(starttime)
         end if

         !
         ! Calculate the sequence of reflectors, updating only a band within A and B
         !

         !$omp task default(none) shared(generationtime, U,V,A,B,tauU,tauV,B_copy) &
         !$omp& private(starttime,endtime,count_rate,istart, istop, jc, work, ib, jb, nb,ierr)&
         !$omp& firstprivate(lda,ldb,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: A) depend(inout: B) depend(out: U)&
         !$omp& depend(out: V)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         U = zero
         V = zero
         tauU = zero
         tauV = zero

         ! ic indicates the column we are reducing this iteration relative to the block
         do ic = 1, nq

            ! jc indicates the column we are reducing this iteration
            jc = jblock + ic - 1

            nb = min( r, n-jc )
            istart = jblock + 1

            !
            ! Reduce a column in A
            !

            !
            ! Apply the previous left reflectors of this block
            ! to B(jblock:jc+r,jc+r).
            !
            if( jc + r .le. n ) then
               do i = 1,ic-1
                  call dlarf( 'L', r, 1, U( i, i, 1 ), 1,&
                     tauU(i, 1), B( jblock + i, jc+r ), ldb, work )
               end do
            end if

            !
            ! Apply the previous left reflectors of this block
            ! to A(jblock+1:jc+r,jc).
            !
            do i = 1,ic-1
               if( jc .le. n ) then
                  call dlarf( 'L', min(r,n-jblock-i+1), 1, U( i, i, 1 ), 1,&
                     tauU(i, 1), A( jblock + i, jc ), lda, work )
               end if
            end do

            if( nb .lt. 2 ) cycle

            ! Reflector to reduce a column in A
            call dlarfg( nb, A( jc+1, jc ), A( jc+2, jc ), 1, tauU(ic, 1) )
            ! Store reflector in U
            do i = 2,nb
               U( ic + i - 1, ic, 1 ) = A( jc+i, jc )
               A( jc+i, jc ) = zero
            end do
            U( ic, ic, 1 ) = one
            ! Apply reflector to B
            call dlarf( 'L', nb, nb, U( ic, ic, 1 ), 1,&
               tauU(ic, 1), B( jc+1, jc+1 ), ldb, work )

            ! RQ factorization of block in B
            call dlacpy( 'A', nb, nb, B( jc+1, jc+1 ), ldb, B_copy, r )
            call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
            ! Form orthogonal matrix
            call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
            ! Reflector to annihilate top row of orthogonal matrix
            call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tauV(ic, 1) )
            do i = 2,nb
               V(ic + i - 1, ic, 1) = B_copy( 1, i )
            end do
            V( ic, ic, 1 ) = one
            ! Apply reflector to small block in A
            call dlarf( 'R', min(n,jc+r+nb)-jc, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), A( jc+1, jc+1 ), lda, work )
            ! Apply reflector to small block in B
            call dlarf( 'R', nb, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), B( jc+1, jc+1 ), ldb, work )
            do i = 2,nb
               B( jc+i, jc+1 ) = zero
            end do
            call dlarf( 'R', jc-istart+1, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), A( istart, jc+1 ), lda, work )
            call dlarf( 'R', jc-istart+1, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), B( istart, jc+1 ), ldb, work )

            !
            ! Chase down the fill in
            !

            ! jb points to the column in A whose fill in is to be reduced in this iteration
            do jb = jc+1, n, r

               ! ib is the block index of the block we are chasing
               ib = (jb - (jc+1))/r + 2
               nb = min( r, n-jb-r+1 )
               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)

               !
               ! Apply the previous left reflectors of this block
               ! to B(jb+r-ic+1:jb+r+r-1,jb+r+r-1).
               !
               if( jb+r+r-1 .le. n ) then
                  do i = 1,ic-1
                     call dlarf( 'L', r, 1, U( i, i, ib ), 1,&
                        tauU(i, ib), B( jb+r-ic + i, jb+r+r-1 ), ldb, work )
                  end do
               end if

               !
               ! Apply the previous left reflectors of this block
               ! to A(jb+r-ic+1:jb+r+nb-1,jb+1).
               !
               do i = 1,ic-1
                  if( jb+r-ic+i .le. n ) call dlarf( 'L', min( r, n - (jb+r-ic+i) + 1 ), 1, &
                     U( i, i, ib ), 1, tauU(i, ib), A( jb+r-ic+i, jb ), lda, work )
               end do

               if( nb .lt. 2 ) cycle

               ! Reflector to reduce a column in A
               call dlarfg( nb, A( jb+r, jb ), A( jb+r+1, jb ), 1, tauU(ic, ib) )
               ! Store reflector in U
               do i = 1,nb-1
                  U( ic+i, ic, ib ) = A( jb+r+i, jb )
                  A( jb+r+i, jb ) = zero
               end do
               U( ic, ic, ib ) = one
               ! Apply reflector to B
               call dlarf( 'L', nb, nb, U( ic, ic, ib ), 1,&
                  tauU(ic, ib), B( jb+r, jb+r ), ldb, work )

               ! RQ factorization of block in B
               call dlacpy( 'A', nb, nb, B( jb+r, jb+r ), ldb, B_copy, r )
               call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
               ! Form orthogonal matrix
               call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
               ! Reflector to annihilate top row of orthogonal matrix
               call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tauV(ic, ib) )
               do i = 2,nb
                  V(ic + i - 1, ic, ib) = B_copy( 1, i )
               end do
               V( ic, ic, ib ) = one
               ! Apply reflector to small block in A
               call dlarf( 'R', min(n,jb+2*r+nb-1)-jb-r+1, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), A( jb+r, jb+r ), lda, work )
               ! Apply reflector to small block in B
               call dlarf( 'R', nb, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), B( jb+r, jb+r ), ldb, work )
               do i = 1,nb-1
                  B( jb+r+i, jb+r ) = zero
               end do
               call dlarf( 'R', jb+r-istart, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), A( istart, jb+r ), lda, work )
               call dlarf( 'R', jb+r-istart, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), B( istart, jb+r ), ldb, work )

            end do

         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            generationtime = generationtime + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task !end of banded update

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            generationtime2 = generationtime2 + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

         !
         ! Lookahead update of B
         !

         !$omp task default(none) shared(V,B,tauV,U,tauU,lookaheadtimeB)&
         !$omp& private(c1,c2,r1,r2,i,starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: B(2,2)) depend(inout: B) depend(in: V) depend(in: U)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on B from the right
         !
         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,V,tauV)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
                  r1 = max( i, istart )
                  r2 = min( i+blocksize-1, istart2-1 )

                  if( r2 .ge. r1 ) then
                     call dlarf( 'R', r2-r1+1, r, V( ic, ic, ib ), 1,&
                        tauV(ic, ib), B( r1, jb+ic ), ldb, work )
                  end if
               end do

            end do
         end do

         !
         ! Lookahead on B from the left
         !
         !$omp taskloop default(none) private(istop,istop2,c1,c2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,U,tauU)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               istop = min( n, jblock + (ib-1)*r  + r + nq - 1 )
               istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )
               c2 = min( istop2, i + blocksize-1 )

               if(c2 .ge. c1) then
                  do ic = 1,nq
                     call dlarf( 'L', min(r,n-jb-ic+1), c2-c1+1, U( ic, ic, ib ), 1,&
                        tauU(ic, ib), B( jb+ic, c1 ), ldb, work )
                  end do
               end if
            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeB = lookaheadtimeB + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         !
         ! Lookahead update of A
         !

         !$omp task default(none) shared(V,A,tauV,U,tauU,lookaheadtimeA)&
         !$omp& private(c1,c2,r1,r2,i,starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: A(2,2)) depend(inout: A) depend(in: V) depend(in: U)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on A from the right
         !
         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,V,tauV)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
                  r1 = max( i, istart )
                  r2 = min( i+blocksize-1, istart2-1 )

                  if( r2 .ge. r1 ) then
                     call dlarf( 'R', r2-r1+1, r, V( ic, ic, ib ), 1,&
                        tauV(ic, ib), A( r1, jb+ic ), lda, work )
                  end if
               end do

            end do

         end do


         !
         ! Lookahead on A from the left
         !
         !$omp taskloop default(none) private(istop,istop2,c1,c2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,U,tauU)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1

               if( ib .eq. 1 ) then
                  jb = jblock
                  istop = min( n, jb  + nq - 1 )
                  istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
                  c1 = max( istop+1, i )
                  c2 = min( istop2, i + blocksize-1 )
                  if(c2 .ge. c1) then
                     do ic = 1,nq
                        call dlarf( 'L', min(r,n-jb-ic+1), c2-c1+1, U( ic, ic, ib ), 1,&
                           tauU(ic, ib), A( jb+ic, c1 ), lda, work )
                     end do
                  end if
               else
                  jb = jblock + 1 + (ib-2)*r
                  istop = min( n, jb  + nq - 1 )
                  istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
                  c1 = max( istop+1, i )
                  c2 = min( istop2, i + blocksize-1 )
                  if(c2 .ge. c1) then
                     do ic = 1,nq
                        call dlarf( 'L', min(r,n-jb-ic-r+2), c2-c1+1, U( ic, ic, ib ), 1,&
                           tauU(ic, ib), A( jb+r+ic-1, c1 ), lda, work )
                     end do
                  end if
               end if

            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeA = lookaheadtimeA + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            lookaheadtime = lookaheadtime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if


         !
         ! Copy V, U, TV and TU to a different location so we can already
         ! start calculating the next reflectors while we update the rest of the matrices
         !
         !$omp task default(none) shared(U,V,U2,V2,tauU,tauU2,tauV,tauV2)&
         !$omp& firstprivate(r,nq,nblocks)&
         !$omp& depend(in: U) depend(in: V)&
         !$omp& depend(out: U2) depend(out: V2)&
         !$omp& priority(100)
         call dcopy( (r+nq)*nq*nblocks, U, 1, U2, 1 )
         call dcopy( (r+nq)*nq*nblocks, V, 1, V2, 1 )
         call dcopy( nq*nblocks, tauV, 1, tauV2, 1 )
         call dcopy( nq*nblocks, tauU, 1, tauU2, 1 )
         !$omp end task

         !
         ! Update B
         !

         !$omp task default(none) shared(V2,B,tauV2,U2,tauU2,updatetimeB)&
         !$omp& private(r1,r2,c1,c2,starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: B(2,2)) depend(in: V2) depend(in: U2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !
         ! Update B from the right
         !

         !$omp taskloop default(none) private(istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,V2,tauV2)
         do i = 1,n,blocksize
            r1 = i

            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  r2 = min( i+blocksize-1, istart2-1 )

                  if(r2 .ge. r1) then
                     call dlarf( 'R', r2-r1+1, min(r,n-jb-ic+1), V2( ic, ic, ib ), 1,&
                        tauV2(ic, ib), B( r1, jb+ic ), ldb, work )
                  end if
               end do
            end do

         end do

         !$omp taskwait


         !
         ! Update B from the left
         !

         !$omp taskloop default(none) private(istop,c1,c2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,U2,tauU2)
         do i = 1,n,blocksize
            c2 = min( n, i + blocksize-1 )

            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )

               if( c2 .ge. c1 ) then

                  do ic = 1,nq
                     call dlarf( 'L', min(r,n-jb-ic+1), c2-c1+1, U2( ic, ic, ib ), 1,&
                        tauU2(ic, ib), B( jb+ic, c1 ), ldb, work )
                  end do
               end if
            end do

         end do

         !$omp taskwait

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeB = updatetimeB + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task


         !
         ! Update A
         !

         !$omp task default(none) shared(V2,A,tauV2,U2,tauU2,updatetimeA)&
         !$omp& private(r1,r2,c1,c2,starttime,endtime, count_rate, work, ib, jb, nb, istart, istop, istart2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: A(2,2)) depend(in: V2) depend(in: U2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)
         !
         ! Update A from the right
         !

         !$omp taskloop default(none) private(istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,V2,tauV2)
         do i = 1,n,blocksize
            r1 = i

            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-ic)*r)
                  r2 = min( i+blocksize-1, istart2-1 )

                  if( r2 .ge. r1 ) then
                     call dlarf( 'R', r2-r1+1, min(r,n-jb-ic+1), V2( ic, ic, ib ), 1,&
                        tauV2(ic, ib), A( r1, jb+ic ), lda, work )
                  end if
               end do
            end do
         end do

         !
         ! Update A from the left
         !

         !$omp taskloop default(none) private(istop,c1,c2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,U2,tauU2)
         do i = 1,n,blocksize
            c2 = min( n, i + blocksize-1 )

            do ib = nblocks,1,-1


               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )

               if( c2 .ge. c1 ) then

                  if( ib .eq. 1 ) then
                     jb = jblock

                     do ic = 1,nq
                        call dlarf( 'L', min(r,n-jb-ic+1), c2-c1+1, U2( ic, ic, ib ), 1,&
                           tauU2(ic, ib), A( jb+ic, c1 ), lda, work )
                     end do
                  else
                     jb = jblock + 1 + (ib-2)*r

                     do ic = 1,nq
                        call dlarf( 'L', min(r,n-jb-ic-r+2), c2-c1+1, U2( ic, ic, ib ), 1,&
                           tauU2(ic, ib), A( jb+r+ic-1, c1 ), lda, work )
                     end do
                  end if

               end if

            end do

         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeA = updatetimeA + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task


         !
         ! Update Q
         !
         !$omp task default(none) shared(U2,Q,tauU2,updatetimeQ)&
         !$omp& private(r1,r2,i, work, ib, jb, nb, starttime, endtime, count_rate)&
         !$omp& firstprivate(ldq,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: Q) depend(in: U2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none) private(r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldq,nblocks,blocksize) &
         !$omp& shared(Q,U2,tauU2)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  call dlarf( 'R', r2-r1+1, min(r,n-jb-ic+1), U2( ic, ic, ib ), 1,&
                     tauU2(ic, ib), Q( r1, jb+ic ), ldq, work )
               end do

            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeQ = updatetimeQ + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task

         !
         ! Update Z
         !
         !$omp task default(none) shared(V2,tauV2,Z,updatetimeZ)&
         !$omp& private(r1,r2,i, work, ib, jb, nb, starttime, endtime, count_rate)&
         !$omp& firstprivate(ldz,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: Z) depend(in: V2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none) private(r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldz,nblocks,blocksize) &
         !$omp& shared(Z,V2,tauV2)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  call dlarf( 'R', r2-r1+1, min(r,n-jb-ic+1), V2( ic, ic, ib ), 1,&
                     tauV2(ic, ib), Z( r1, jb+ic ), ldz, work )
               end do
            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeZ = updatetimeZ + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            applicationtime = applicationtime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

      end do

      !$omp taskwait
      !$omp end single nowait
      !$omp end parallel

      if(time_tasks) then
         write(*,'("      gen: ",F16.6," ;lookahead A: ", F16.6," ;lookahead B: ", F16.6,&
         &" ;update A: ", F16.6," ;update B: ", F16.6," ;update Q: ", F16.6," ;update Z: ", F16.6)')&
            generationtime, lookaheadtimeA, lookaheadtimeB, updatetimeA, updatetimeB, updatetimeQ, updatetimeZ
      end if
      if(time_components) then
         write(*,'("      gen: ",F16.6," ;lookahead: ", F16.6," ;application: ", F16.6)')&
            generationtime2, lookaheadtime, applicationtime
      end if

      deallocate( U, V, tauU, tauV, work, B_copy, U2, V2 )

   end subroutine

   subroutine rht2ht_givens( n, A, lda, B, ldb, Q, ldq, Z, ldz, r, nq )
      implicit none
      ! arguments
      integer, intent(in) :: n, lda, ldb, ldq, ldz, r, nq
      double precision, intent(inout) :: A(lda,*), B(ldb,*), Z(ldz,*),Q(ldq,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, parameter :: min_blocksize = 256
      ! Number of accumulated rotation matrices to combine into a larger matrix to improve performance
      integer, parameter :: mb = 1
      logical, parameter :: parallelize_dgemm = .false.

      ! local variables
      integer :: i,ig,jc, jb, ierr, nb, nblocks, jblock, ic, ib, istart, istop, istart2,&
         istop2,ibm,ibm2,nthreads,blocksize
      double precision, allocatable :: work(:), U(:,:,:), V(:,:,:), &
         CU(:,:,:), SU(:,:,:), CV(:,:,:), SV(:,:,:),&
         U2(:,:,:), V2(:,:,:), CV2(:,:,:), SV2(:,:,:),U3(:,:,:),V3(:,:,:)
      double precision :: temp

      ! temporary timing variables
      logical, parameter :: time_tasks = .true.
      integer(8) :: starttime,endtime,count_rate
      double precision :: generationtime, applicationtime, lookaheadtimeA, lookaheadtimeB,&
         updatetimeB

      ! ---------------
      ! executable code
      ! ---------------

      nblocks = (n-1)/r
      if( mod( (n-1),r ) .ne. 0 ) nblocks = nblocks + 1


      allocate(work( (mb*r+nq)*n ))
      allocate( U( r+nq,  r+nq, nblocks ) )
      allocate( V( r+nq,  r+nq, nblocks ) )
      allocate( CU( r+nq, nq, nblocks ) )
      allocate( SU( r+nq, nq, nblocks ) )
      allocate( CV( r+nq, nq, nblocks ) )
      allocate( SV( r+nq, nq, nblocks ) )
      allocate( U2( r+nq,  r+nq, nblocks ) )
      allocate( V2( r+nq,  r+nq, nblocks ) )
      allocate( CV2( r+nq, nq, nblocks ) )
      allocate( SV2( r+nq, nq, nblocks ) )
      allocate( U3( mb*r+nq,  mb*r+nq, nblocks ) )
      allocate( V3( mb*r+nq,  mb*r+nq, nblocks ) )

      !$omp parallel
      !$omp single

      generationtime = 0.0d0
      applicationtime = 0.0d0
      lookaheadtimeA = 0.0d0
      lookaheadtimeB = 0.0d0

      nthreads = omp_get_max_threads()
      blocksize = n/nthreads
      if( mod( blocksize, 16 ) .ne. 0 ) blocksize = blocksize + (16-mod( blocksize, 16 ))
      blocksize = max( min_blocksize, blocksize )
      blocksize = min( n, blocksize )

      ! jblock indicates the first column of the block we are reducing
      do jblock = 1,n-2,nq
         ! do jblock = 1,1,nq

         nblocks = (n-jblock)/r
         if( mod( (n-jblock),r ) .ne. 0 ) nblocks = nblocks + 1


         !
         ! Calculate the sequence of rotations, updating only a band within A and B
         !

         !$omp task default(none) shared(generationtime, CU,SU,CV,SV,A,B) &
         !$omp& private(starttime,endtime,count_rate,istart, istop, jc, ib, jb, nb,ierr,temp)&
         !$omp& firstprivate(lda,ldb,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: A) depend(inout: B)&
         !$omp& depend(out: CU) depend(out: SU) depend(out: CV) depend(out: SV)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         CU = one
         CV = one
         SU = zero
         SV = zero

         ! ic indicates the column we are reducing this iteration relative to the block
         do ic = 1, nq

            ! jc indicates the column we are reducing this iteration
            jc = jblock + ic - 1

            nb = min( r, n-jc )
            istart = jblock + 1

            !
            ! Reduce a column in A
            !

            !
            ! Apply the previous left reflectors of this block
            ! to B(jblock:jc+r,jc+r).
            !
            if( jc + r .le. n ) then
               do i = 1,ic-1
                  do ig = r-1,1,-1
                     call drot( 1, B( jblock+i-1+ig, jc+r ),ldb,&
                        B( jblock+i+ig,jc+r),ldb,&
                        CU(i+ig-1,i,1), SU(i+ig-1,i,1) )
                  end do
               end do
            end if

            !
            ! Apply the previous left reflectors of this block
            ! to A(jblock+1:jc+r,jc).
            !
            do i = 1,ic-1
               if( jc .le. n ) then
                  do ig = min(r-1,n - (jblock+i)),1,-1
                     call drot( 1, A( jblock+i-1+ig, jc ),lda,&
                        A( jblock+i+ig,jc),lda,&
                        CU(i+ig-1,i,1), SU(i+ig-1,i,1) )
                  end do
               end if
            end do

            if( nb .lt. 2 ) cycle

            ! Rotations to reduce a column in A
            do ig = nb-1,1,-1
               call dlartg( A( jc+ig, jc ), A( jc+ig+1, jc ),&
                  CU(ic+ig-1,ic,1), SU(ic+ig-1,ic,1), temp )
               A( jc+ig, jc ) = temp
               A( jc+ig+1, jc ) = zero
            end do
            ! Apply rotations to B
            do ig = nb-1,1,-1
               call drot( nb, B( jc+ig, jc+1 ),ldb,B( jc+ig+1,jc+1),ldb,&
                  CU(ic+ig-1,ic,1), SU(ic+ig-1,ic,1) )
            end do

            ! Rotations to reduce the fill in in B
            do ig = nb-1,1,-1
               call dlartg( B( jc+ig+1, jc+ig+1 ), B( jc+ig+1, jc+ig ),&
                  CV(ic+ig-1,ic,1), SV(ic+ig-1,ic,1), temp )
               B( jc+ig+1, jc+ig+1 ) = temp
               B( jc+ig+1, jc+ig ) = zero

               call drot( ig, B( jc+1, jc+ig+1 ),1,B( jc+1,jc+ig),1,&
                  CV(ic+ig-1,ic,1), SV(ic+ig-1,ic,1) )
            end do
            ! Apply rotations to small block in A
            do ig = nb-1,1,-1
               call drot( min(n,jc+r+nb)-jc, A( jc+1, jc+ig+1 ),1,A( jc+1,jc+ig),1,&
                  CV(ic+ig-1,ic,1), SV(ic+ig-1,ic,1) )
            end do
            ! Apply rotations to block in B
            !$omp task
            do ig = nb-1,1,-1
               call drot( jc-istart+1, B( istart, jc+ig+1 ),1,B( istart,jc+ig),1,&
                  CV(ic+ig-1,ic,1), SV(ic+ig-1,ic,1) )
            end do
            ! Apply rotations to block in A
            do ig = nb-1,1,-1
               call drot( jc-istart+1, A( istart, jc+ig+1 ),1,A( istart,jc+ig),1,&
                  CV(ic+ig-1,ic,1), SV(ic+ig-1,ic,1) )
            end do
            !$omp end task

            !
            ! Chase down the fill in
            !

            ! jb points to the column in A whose fill in is to be reduced in this iteration
            do jb = jc+1, n, r
               ! do jb = jc+1, jc+r+1, r

               ! ib is the block index of the block we are chasing
               ib = (jb - (jc+1))/r + 2
               nb = min( r, n-jb-r+1 )
               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)

               !
               ! Apply the previous left reflectors of this block
               ! to B(jb+r-ic+1:jb+r+r-1,jb+r+r-1).
               !
               if( jb+r+r-1 .le. n ) then
                  do i = 1,ic-1
                     do ig = r-1,1,-1
                        call drot( 1, B( jb+r-ic+i+ig-1,jb+r+r-1 ),ldb,&
                           B( jb+r-ic+i+ig,jb+r+r-1),ldb,&
                           CU(i+ig-1,i,ib), SU(i+ig-1,i,ib) )
                     end do
                  end do
               end if

               !
               ! Apply the previous left reflectors of this block
               ! to A(jb+r-ic+1:jb+r+nb-1,jb+1).
               !
               do i = 1,ic-1
                  if( jb+r-ic+i .le. n ) then
                     do ig = min(r,n - (jb+r-ic+i) + 1)-1,1,-1
                        call drot( 1, A( jb+r-ic+i+ig-1,jb ),lda,&
                           A( jb+r-ic+i+ig,jb),lda,&
                           CU(i+ig-1,i,ib), SU(i+ig-1,i,ib) )
                     end do
                  end if
               end do

               if(nb .lt. 2) cycle

               ! Rotations to reduce a column in A
               do ig = nb-1,1,-1
                  call dlartg( A( jb+r-1+ig, jb ), A( jb+r+ig, jb ),&
                     CU(ic+ig-1,ic,ib), SU(ic+ig-1,ic,ib), temp )
                  A( jb+ig+r-1, jb ) = temp
                  A( jb+ig+r, jb ) = zero
               end do
               ! Apply rotations to B
               do ig = nb-1,1,-1
                  call drot( nb, B( jb+r-1+ig, jb+r ),ldb,B( jb+ig+r,jb+r),ldb,&
                     CU(ic+ig-1,ic,ib), SU(ic+ig-1,ic,ib) )
               end do

               ! Rotations to reduce the fill in in B
               do ig = nb-1,1,-1
                  call dlartg( B( jb+r+ig, jb+r+ig ), B( jb+r+ig, jb+r+ig-1 ),&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib), temp )
                  B( jb+r+ig, jb+r+ig ) = temp
                  B( jb+r+ig, jb+r+ig-1 ) = zero

                  call drot( ig, B( jb+r, jb+r+ig ),1,B( jb+r,jb+r+ig-1),1,&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
               end do
               ! Apply rotations to small block in A
               do ig = nb-1,1,-1
                  call drot( min(n,jb+2*r+nb-1)-(jb+r)+1, A( jb+r, jb+r+ig ),1,A( jb+r,jb+r+ig-1),1,&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
               end do
               ! Apply rotations to block in B
               !$omp task
               do ig = nb-1,1,-1
                  call drot( jb+r-istart, B( istart, jb+r+ig ),1,B( istart,jb+r+ig-1),1,&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
               end do
               ! Apply rotations to block in A
               do ig = nb-1,1,-1
                  call drot( jb+r-istart, A( istart, jb+r+ig ),1,A( istart,jb+r+ig-1),1,&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
               end do
               !$omp end task

            end do
            !$omp taskwait

         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            generationtime = generationtime + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task


         !
         ! Accumulate left rotations
         !
         !$omp task default(none) shared(CU,SU,U)&
         !$omp& private(ib, ic, ig)&
         !$omp& firstprivate(n,r,nq,nblocks)&
         !$omp& depend(in: CU) depend(in: SU) depend(out: U)&
         !$omp& priority(100)
         do ib = 1,nblocks
            call dlaset( 'A', r+nq,r+nq,zero,one,U(1,1,ib),r+nq)
            do ic = 1, nq
               do ig = r-1,1,-1
                  call drot( r+nq, U(1,ic-1+ig,ib),1,U(1,ic+ig,ib),1,&
                     CU(ic+ig-1,ic,ib), SU(ic+ig-1,ic,ib) )
               end do
            end do
         end do
         !$omp end task

         !
         ! Accumulate right rotations
         !
         !$omp task default(none) shared(CV,SV,V)&
         !$omp& private(ib, ic, ig)&
         !$omp& firstprivate(n,r,nq,nblocks)&
         !$omp& depend(in: CV) depend(in: SV) depend(out: V)&
         !$omp& priority(100)
         do ib = 1,nblocks
            call dlaset( 'A', r+nq,r+nq,zero,one,V(1,1,ib),r+nq)
            do ic = 1, nq
               do ig = r-1,1,-1
                  call drot( r+nq, V(1,ic+ig,ib),1,V(1,ic+ig-1,ib),1,&
                     CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
               end do
            end do
         end do
         !$omp end task

         !
         ! Lookahead updates of B
         !
         !$omp task default(none) shared(SV,CV,U,V,B,lookaheadtimeB)&
         !$omp& private(starttime, endtime, count_rate, ib, jb, nb, ig, work, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: B(2,2)) depend(inout: B) depend(in: V) depend(in: U) depend(in: CV) depend(in: SV)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on B from the right
         !
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )

            do ic = 1,nq
               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
               istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
               if(istart2 > istart) then
                  do ig = r-1,1,-1
                     call drot( istart2-istart, B( istart, jb+ic+ig ),1,B( istart,jb+ic+ig-1),1,&
                        CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
                  end do
               end if
            end do

         end do

         !
         ! Lookahead on B from the left
         !
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )

            istop = min( n, jblock + (ib-1)*r  + r + nq - 1 )
            istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
            if( istop2 > istop ) then
               call dgemm( 'T', 'N', nb, istop2-istop, nb, one, U(1,1,ib), r+nq, B(jb+1,istop+1), ldb, zero, work, nb )
               call dlacpy( 'A',nb,istop2-istop,work,nb,B(jb+1,istop+1),ldb )
            end if
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeB = lookaheadtimeB + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task

         !
         ! Lookahead updates of A
         !
         !$omp task default(none) shared(SV,CV,U,V,A,lookaheadtimeA)&
         !$omp& private(starttime,endtime,count_rate,ib, jb, nb, ig, work, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: A(2,2)) depend(inout: A) depend(in: V) depend(in: U) depend(in: CV) depend(in: SV)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on A from the right
         !
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )


            do ic = 1,nq
               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
               istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
               if(istart2 > istart) then
                  do ig = min(r-1,n-jb-ic),1,-1
                     call drot( istart2-istart, A( istart, jb+ic+ig ),1,A( istart,jb+ic+ig-1),1,&
                        CV(ic+ig-1,ic,ib), SV(ic+ig-1,ic,ib) )
                  end do
               end if
            end do

         end do

         !
         ! Lookahead on A from the left
         !
         do ib = nblocks,1,-1

            if( ib .eq. 1 ) then
               jb = jblock
               nb = min( n-jb, r+nq )
               istop = min( n, jb  + nq - 1 )
               istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               if( istop2 > istop ) then
                  call dgemm( 'T', 'N', nb, istop2-istop, nb, one, U(1,1,ib), r+nq, A(jb+1,istop+1), lda, zero, work, nb )
                  call dlacpy( 'A',nb,istop2-istop,work,nb,A(jb+1,istop+1),lda )
               end if
            else
               jb = jblock + 1 + (ib-2)*r
               nb = min( n-jb-r+1, r+nq )
               istop = min( n, jb  + nq - 1 )
               istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               if( istop2 > istop ) then
                  call dgemm( 'T', 'N', nb, istop2-istop, nb, one, U(1,1,ib), r+nq, A(jb+r,istop+1), lda, zero, work, nb )
                  call dlacpy( 'A',nb,istop2-istop,work,nb,A(jb+r,istop+1),lda )
               end if
            end if

         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeA = lookaheadtimeA + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         !
         ! Copy U, V, CV and SV to a different location so we can already
         ! start calculating the next rotations while we update the rest of the matrices
         !
         !$omp task default(none) shared(U,V,CV,SV,U2,V2,CV2,SV2)&
         !$omp& firstprivate(r,nq,nblocks)&
         !$omp& depend(in: U) depend(in: V) depend(in: CV) depend(in: SV)&
         !$omp& depend(out: U2) depend(out: V2) depend(out: CV2) depend(out: SV2)&
         !$omp& priority(100)
         call dcopy( (r+nq)*(r+nq)*nblocks, U, 1, U2, 1 )
         call dcopy( (r+nq)*(r+nq)*nblocks, V, 1, V2, 1 )
         call dcopy( (r+nq)*nq*nblocks, CV, 1, CV2, 1 )
         call dcopy( (r+nq)*nq*nblocks, SV, 1, SV2, 1 )
         !$omp end task

         !
         ! Update the rest of B
         !
         !$omp task default(none) shared(SV2,CV2,U2,V2,B,V3,updatetimeB)&
         !$omp& private(starttime,endtime,count_rate,ibm, ibm2, i, ib, jb, nb, ig, work, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: B(2,2)) depend(in: V2) depend(in: U2) depend(in: CV2) depend(in: SV2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !
         ! Update B from the right
         !
         if( mb .gt. 1 ) then

            ibm2 = nblocks
            ibm = max(nblocks-mb+1,1)
            do while(ibm .ge. 1)

               ! Perform unblocked updates from the right up to the edge of the first block in the multiblock structure
               istart = max(jblock + 1,jblock + 1 + (ibm-1)*r - (2*nq-1)*r)
               do ib = ibm2,ibm,-1
                  jb = jblock + (ib-1)*r
                  nb = min( n-jb, r+nq )


                  do ic = 1,nq
                     istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                     if(istart2 > istart) then
                        do ig = r-1,1,-1
                           call drot( istart2-istart, B( istart, jb+ic+ig ),1,B( istart,jb+ic+ig-1),1,&
                              CV2(ic+ig-1,ic,ib), SV2(ic+ig-1,ic,ib) )
                        end do
                     end if
                  end do

               end do

               do ib = ibm2,ibm,-1
                  jb = jblock + (ib-1)*r
                  nb = min( n-jb, r+nq )

                  call dgemm( 'N', 'N', istart-1, nb, nb, one, B(1,jb+1), ldb, V2(1,1,ib), r+nq, zero, work, istart-1 )
                  call dlacpy( 'A', istart-1, nb, work, istart-1, B(1,jb+1),ldb )

               end do
               if(ibm .eq. 1) exit
               ibm2 = ibm-1
               ibm = max(ibm-mb,1)
            end do

         else

            !
            ! Update B from the right
            !
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )

               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-1)*r)
               do ic = 1,nq
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  if(istart2 > istart) then
                     do ig = r-1,1,-1
                        call drot( istart2-istart, B( istart, jb+ic+ig ),1,B( istart,jb+ic+ig-1),1,&
                           CV2(ic+ig-1,ic,ib), SV2(ic+ig-1,ic,ib) )
                     end do
                  end if
               end do

               call dgemm( 'N', 'N', istart-1, nb, nb, one, B(1,jb+1), ldb, V2(1,1,ib), r+nq, zero, work, istart-1 )
               call dlacpy( 'A', istart-1, nb, work, istart-1, B(1,jb+1),ldb )
            end do

         end if

         !
         ! Update B from the left
         !
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )

            istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
            if( n > istop ) then

               call dgemm( 'T', 'N', nb, n-istop, nb, one, U2(1,1,ib), r+nq, B(jb+1,istop+1), ldb, zero, work, nb )
               call dlacpy( 'A',nb,n-istop,work,nb,B(jb+1,istop+1),ldb )
            end if
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeB = updatetimeB + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task


         !
         ! Update the rest of A
         !
         !$omp task default(none) shared(SV2,CV2,U2,V2,A)&
         !$omp& private(ib, jb, nb, ig, work, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: A(2,2)) depend(in: V2) depend(in: U2) depend(in: CV2) depend(in: SV2)&
         !$omp& priority(1)

         !
         ! Update A from the right
         !
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )

            ! Perform the right updates up to the edge of the block
            istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-1)*r)
            do ic = 1,nq
               istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-ic)*r)
               do ig = min(r-1,n-jb-ic),1,-1
                  call drot( istart2-istart, A( istart, jb+ic+ig ),1,A( istart,jb+ic+ig-1),1,&
                     CV2(ic+ig-1,ic,ib), SV2(ic+ig-1,ic,ib) )
               end do
            end do

            call dgemm( 'N', 'N', istart-1, nb, nb, one, A(1,jb+1), lda, V2(1,1,ib), r+nq, zero, work, istart-1 )
            call dlacpy( 'A', istart-1, nb, work, istart-1, A(1,jb+1),lda )
         end do

         !
         ! Update A from the left
         !
         do ib = nblocks,1,-1

            if( ib .eq. 1 ) then
               jb = jblock
               nb = min( n-jb, r+nq )
               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               if( n > istop ) then
                  call dgemm( 'T', 'N', nb, n-istop, nb, one, U2(1,1,ib), r+nq, A(jb+1,istop+1), lda, zero, work, nb )
                  call dlacpy( 'A',nb,n-istop,work,nb,A(jb+1,istop+1),lda )
               end if
            else
               jb = jblock + 1 + (ib-2)*r
               nb = min( n-jb-r+1, r+nq )
               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               if( n > istop ) then
                  call dgemm( 'T', 'N', nb, n-istop, nb, one, U2(1,1,ib), r+nq, A(jb+r,istop+1), lda, zero, work, nb )
                  call dlacpy( 'A',nb,n-istop,work,nb,A(jb+r,istop+1),lda )
               end if
            end if

         end do
         !$omp end task


         !
         ! Update Q
         !
         !$omp task default(none) shared(U2,Q)&
         !$omp& private(i, ib, jb, nb, work)&
         !$omp& firstprivate(ldq,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: Q) depend(in: U2)&
         !$omp& priority(1)

         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )
            call dgemm( 'N', 'N', n, nb, nb, one, Q(1,jb+1), ldq, U2(1,1,ib), r+nq, zero, work, n )
            call dlacpy( 'A', n, nb, work, n, Q(1,jb+1), ldq )
         end do
         !$omp end task

         !
         ! Update Z
         !
         !$omp task default(none) shared(V2,Z)&
         !$omp& private(i, ib, jb, nb, work)&
         !$omp& firstprivate(ldz,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: Z) depend(in: V2)&
         !$omp& priority(1)
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )
            call dgemm( 'N', 'N', n, nb, nb, one, Z(1,jb+1), ldz, V2(1,1,ib), r+nq, zero, work, n )
            call dlacpy( 'A', n, nb, work, n, Z(1,jb+1), ldz )
         end do
         !$omp end task

      end do

      !$omp taskwait
      !$omp end single nowait
      !$omp end parallel


      if(time_tasks) then
         write(*,*) "generate time", generationtime
         write(*,*) "lookahead A time", lookaheadtimeA
         write(*,*) "lookahead B time", lookaheadtimeB
         write(*,*) "update B time", updatetimeB
      end if

      deallocate( U, V, CU, CV, SU, SV, U2, V2, CV2, SV2, work )

   end subroutine

   subroutine rht2ht_householder( n, A, lda, B, ldb, Q, ldq, Z, ldz, r, nq, time_components )
      implicit none
      ! arguments
      integer, intent(in) :: n, lda, ldb, ldq, ldz, r, nq
      double precision, intent(inout) :: A(lda,*), B(ldb,*), Z(ldz,*),Q(ldq,*)
      logical :: time_components

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, parameter :: min_blocksize = 32

      ! local variables
      integer :: r1,r2,c1,c2,i, jc, jb, ierr, nb, nblocks, jblock, ic, ib, istart, istop, istart2,&
         istop2, nthreads, blocksize
      double precision, allocatable :: U(:,:,:), V(:,:,:), TU(:,:,:), TV(:,:,:), &
         work(:), B_copy(:,:), tauU(:,:), tauV(:,:),&
         U2(:,:,:), V2(:,:,:), TU2(:,:,:), TV2(:,:,:),tauU2(:,:), tauV2(:,:)

      ! temporary timing variables
      logical, parameter :: time_tasks = .false.
      integer(8) :: starttime,endtime,count_rate
      double precision :: generationtime, lookaheadtimeA, lookaheadtimeB,&
         updatetimeQ, updatetimeA, updatetimeB, updatetimeZ
      double precision :: generationtime2, applicationtime, lookaheadtime

      ! ---------------
      ! executable code
      ! ---------------


      nthreads = omp_get_max_threads()
      blocksize = n/nthreads
      if( mod( blocksize, 16 ) .ne. 0 ) blocksize = blocksize + (16-mod( blocksize, 16 ))
      blocksize = max( min_blocksize, blocksize )
      blocksize = min( n, blocksize )

      nblocks = (n-1)/r
      if( mod( (n-1),r ) .ne. 0 ) nblocks = nblocks + 1


      allocate(work( (r+nq)*n ))
      allocate(B_copy( r, r ))
      allocate( U( r+nq, nq, nblocks ) )
      allocate( V( r+nq, nq, nblocks ) )
      allocate( TU( nq, nq, nblocks ) )
      allocate( TV( nq, nq, nblocks ) )
      allocate( tauU( nq, nblocks ) )
      allocate( tauV( nq, nblocks ) )
      allocate( U2( r+nq, nq, nblocks ) )
      allocate( V2( r+nq, nq, nblocks ) )
      allocate( TU2( nq, nq, nblocks ) )
      allocate( TV2( nq, nq, nblocks ) )
      allocate( tauU2( nq, nblocks ) )
      allocate( tauV2( nq, nblocks ) )

      !$omp parallel
      !$omp single

      generationtime = 0.0d0
      applicationtime = 0.0d0
      generationtime2 = 0.0d0
      lookaheadtime = 0.0d0
      lookaheadtimeA = 0.0d0
      lookaheadtimeB = 0.0d0
      updatetimeA = 0.0d0
      updatetimeB = 0.0d0
      updatetimeQ = 0.0d0
      updatetimeZ = 0.0d0

      ! jblock indicates the first column of the block we are reducing
      do jblock = 1,n-2,nq
         ! do jblock = 1,1,nq

         nblocks = (n-jblock)/r
         if( mod( (n-jblock),r ) .ne. 0 ) nblocks = nblocks + 1

         if( time_components ) then
            !$omp taskwait
            call system_clock(starttime)
         end if


         !
         ! Calculate the sequence of reflectors, updating only a band within A and B
         !

         !$omp task default(none) shared(generationtime, U,TU,V,TV,A,B,tauU,tauV,B_copy) &
         !$omp& private(starttime,endtime,count_rate,istart, istop, jc, work, ib, jb, nb,ierr)&
         !$omp& firstprivate(lda,ldb,n,r,nq,jblock,nblocks)&
         !$omp& depend(inout: A) depend(inout: B) depend(out: U) depend(out: TU)&
         !$omp& depend(out: V) depend(out: TV)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         U = zero
         V = zero
         TU = zero
         TV = zero
         tauU = zero
         tauV = zero

         ! ic indicates the column we are reducing this iteration relative to the block
         do ic = 1, nq

            ! jc indicates the column we are reducing this iteration
            jc = jblock + ic - 1

            nb = min( r, n-jc )
            istart = jblock + 1

            !
            ! Reduce a column in A
            !

            !
            ! Apply the previous left reflectors of this block
            ! to B(jblock:jc+r,jc+r).
            !
            if( jc + r .le. n ) then
               do i = 1,ic-1
                  call dlarf( 'L', r, 1, U( i, i, 1 ), 1,&
                     tauU(i, 1), B( jblock + i, jc+r ), ldb, work )
               end do
            end if

            !
            ! Apply the previous left reflectors of this block
            ! to A(jblock+1:jc+r,jc).
            !
            do i = 1,ic-1
               if( jc .le. n ) then
                  call dlarf( 'L', min(r,n-jblock-i+1), 1, U( i, i, 1 ), 1,&
                     tauU(i, 1), A( jblock + i, jc ), lda, work )
               end if
            end do

            if( nb .lt. 2 ) cycle

            ! Reflector to reduce a column in A
            call dlarfg( nb, A( jc+1, jc ), A( jc+2, jc ), 1, tauU(ic, 1) )
            ! Store reflector in U
            do i = 2,nb
               U( ic + i - 1, ic, 1 ) = A( jc+i, jc )
               A( jc+i, jc ) = zero
            end do
            U( ic, ic, 1 ) = one
            ! Apply reflector to B
            call dlarf( 'L', nb, nb, U( ic, ic, 1 ), 1,&
               tauU(ic, 1), B( jc+1, jc+1 ), ldb, work )

            ! RQ factorization of block in B
            call dlacpy( 'A', nb, nb, B( jc+1, jc+1 ), ldb, B_copy, r )
            call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
            ! Form orthogonal matrix
            call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
            ! Reflector to annihilate top row of orthogonal matrix
            call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tauV(ic, 1) )
            do i = 2,nb
               V(ic + i - 1, ic, 1) = B_copy( 1, i )
            end do
            V( ic, ic, 1 ) = one
            ! Apply reflector to small block in A
            call dlarf( 'R', min(n,jc+r+nb)-jc, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), A( jc+1, jc+1 ), lda, work )
            ! Apply reflector to small block in B
            call dlarf( 'R', nb, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), B( jc+1, jc+1 ), ldb, work )
            do i = 2,nb
               B( jc+i, jc+1 ) = zero
            end do
            call dlarf( 'R', jc-istart+1, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), A( istart, jc+1 ), lda, work )
            call dlarf( 'R', jc-istart+1, nb, V( ic, ic, 1 ), 1,&
               tauV(ic, 1), B( istart, jc+1 ), ldb, work )

            !
            ! Chase down the fill in
            !

            ! jb points to the column in A whose fill in is to be reduced in this iteration
            do jb = jc+1, n, r

               ! ib is the block index of the block we are chasing
               ib = (jb - (jc+1))/r + 2
               nb = min( r, n-jb-r+1 )
               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)

               !
               ! Apply the previous left reflectors of this block
               ! to B(jb+r-ic+1:jb+r+r-1,jb+r+r-1).
               !
               if( jb+r+r-1 .le. n ) then
                  do i = 1,ic-1
                     call dlarf( 'L', r, 1, U( i, i, ib ), 1,&
                        tauU(i, ib), B( jb+r-ic + i, jb+r+r-1 ), ldb, work )
                  end do
               end if

               !
               ! Apply the previous left reflectors of this block
               ! to A(jb+r-ic+1:jb+r+nb-1,jb+1).
               !
               do i = 1,ic-1
                  if( jb+r-ic+i .le. n ) call dlarf( 'L', min( r, n - (jb+r-ic+i) + 1 ), 1, &
                     U( i, i, ib ), 1, tauU(i, ib), A( jb+r-ic+i, jb ), lda, work )
               end do

               if( nb .lt. 2 ) cycle

               ! Reflector to reduce a column in A
               call dlarfg( nb, A( jb+r, jb ), A( jb+r+1, jb ), 1, tauU(ic, ib) )
               ! Store reflector in U
               do i = 1,nb-1
                  U( ic+i, ic, ib ) = A( jb+r+i, jb )
                  A( jb+r+i, jb ) = zero
               end do
               U( ic, ic, ib ) = one
               ! Apply reflector to B
               call dlarf( 'L', nb, nb, U( ic, ic, ib ), 1,&
                  tauU(ic, ib), B( jb+r, jb+r ), ldb, work )

               ! RQ factorization of block in B
               call dlacpy( 'A', nb, nb, B( jb+r, jb+r ), ldb, B_copy, r )
               call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
               ! Form orthogonal matrix
               call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
               ! Reflector to annihilate top row of orthogonal matrix
               call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tauV(ic, ib) )
               do i = 2,nb
                  V(ic + i - 1, ic, ib) = B_copy( 1, i )
               end do
               V( ic, ic, ib ) = one
               ! Apply reflector to small block in A
               call dlarf( 'R', min(n,jb+2*r+nb-1)-jb-r+1, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), A( jb+r, jb+r ), lda, work )
               ! Apply reflector to small block in B
               call dlarf( 'R', nb, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), B( jb+r, jb+r ), ldb, work )
               do i = 1,nb-1
                  B( jb+r+i, jb+r ) = zero
               end do
               call dlarf( 'R', jb+r-istart, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), A( istart, jb+r ), lda, work )
               call dlarf( 'R', jb+r-istart, nb, V( ic, ic, ib ), 1,&
                  tauV(ic, ib), B( istart, jb+r ), ldb, work )

            end do

         end do

         ! Form TU
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )
            call dlarft( 'F', 'C', nb, nq, U(1,1,ib), r+nq, tauU( 1, ib ), TU(1,1,ib),nq )
         end do

         ! Form TV
         do ib = nblocks,1,-1
            jb = jblock + (ib-1)*r
            nb = min( n-jb, r+nq )
            call dlarft( 'F', 'C', nb, nq, V(1,1,ib), r+nq, tauV( 1, ib ), TV(1,1,ib),nq )
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            generationtime = generationtime + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task !end of banded update

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            generationtime2 = generationtime2 + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

         !
         ! Lookahead update of B
         !

         !$omp task default(none) shared(V,TV,B,tauV,U,TU,lookaheadtimeB,tauU)&
         !$omp& private(starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: B(2,2)) depend(inout: B) depend(in: V) depend(in: TV) depend(in: U) depend(in: TU)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on B from the right
         !
         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,V,tauV)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
                  r1 = max( i, istart )
                  r2 = min( i+blocksize-1, istart2-1 )

                  if( r2 .ge. r1 ) then
                     call dlarf( 'R', r2-r1+1, r, V( ic, ic, ib ), 1,&
                        tauV(ic, ib), B( r1, jb+ic ), ldb, work )
                  end if
               end do

            end do
         end do

         !
         ! Lookahead on B from the left
         !
         !$omp taskloop default(none) private(istop,istop2,c1,c2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,U,TU)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )

               istop = min( n, jblock + (ib-1)*r  + r + nq - 1 )
               istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )
               c2 = min( istop2, i + blocksize-1 )

               if(c2 .ge. c1) then
                  call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, min(nq,nb), U(1,1,ib), r+nq, TU(1,1,ib),nq, &
                     B(jb+1, c1), ldb, work, c2-c1+1 )
               end if
            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeB = lookaheadtimeB + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         !
         ! Lookahead update of A
         !

         !$omp task default(none) shared(V,TV,A,tauV,U,TU,lookaheadtimeA,tauU)&
         !$omp& private(starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, istop2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: A(2,2)) depend(inout: A) depend(in: V) depend(in: TV) depend(in: U) depend(in: TU)&
         !$omp& priority(100)

         if(time_tasks) call system_clock(starttime)

         !
         ! Lookahead on A from the right
         !
         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,V,tauV)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r

               do ic = 1,nq
                  istart = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                  istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (nq-ic)*r)
                  r1 = max( i, istart )
                  r2 = min( i+blocksize-1, istart2-1 )

                  if( r2 .ge. r1 ) then
                     call dlarf( 'R', r2-r1+1, r, V( ic, ic, ib ), 1,&
                        tauV(ic, ib), A( r1, jb+ic ), lda, work )
                  end if
               end do

            end do

         end do


         !
         ! Lookahead on A from the left
         !
         !$omp taskloop default(none) private(istop,istop2,c1,c2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,U,tauU,TU)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1

               if( ib .eq. 1 ) then
                  jb = jblock
                  nb = min( n-jb, r+nq )
                  istop = min( n, jb  + nq - 1 )
                  istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
                  c1 = max( istop+1, i )
                  c2 = min( istop2, i + blocksize-1 )
                  if(c2 .ge. c1) then
                     call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, nq, U(1,1,ib), r+nq, TU(1,1,ib),nq, &
                        A(jb+1, c1), lda, work, c2-c1+1 )
                  end if
               else
                  jb = jblock + 1 + (ib-2)*r
                  nb = min( n-jb-r+1, r+nq )
                  istop = min( n, jb  + nq - 1 )
                  istop2 = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
                  c1 = max( istop+1, i )
                  c2 = min( istop2, i + blocksize-1 )
                  if(c2 .ge. c1) then
                     call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, nq, U(1,1,ib), r+nq, TU(1,1,ib),nq, &
                        A(jb+r, c1), lda, work, c2-c1+1 )
                  end if
               end if

            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            lookaheadtimeA = lookaheadtimeA + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            lookaheadtime = lookaheadtime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if


         !
         ! Copy V, U, TV and TU to a different location so we can already
         ! start calculating the next reflectors while we update the rest of the matrices
         !
         !$omp task default(none) shared(U,TU,V,TV,U2,TU2,V2,TV2,tauU,tauU2,tauV,tauV2)&
         !$omp& firstprivate(r,nq,nblocks)&
         !$omp& depend(in: U) depend(in: V) depend(in: TU) depend(in: TV)&
         !$omp& depend(out: U2) depend(out: V2) depend(out: TU2) depend(out: TV2)&
         !$omp& priority(100)
         call dcopy( (r+nq)*nq*nblocks, U, 1, U2, 1 )
         call dcopy( (r+nq)*nq*nblocks, V, 1, V2, 1 )
         call dcopy( nq*nq*nblocks, TU, 1, TU2, 1 )
         call dcopy( nq*nq*nblocks, TV, 1, TV2, 1 )
         call dcopy( nq*nblocks, tauV, 1, tauV2, 1 )
         call dcopy( nq*nblocks, tauU, 1, tauU2, 1 )
         !$omp end task

         !
         ! Update B
         !

         !$omp task default(none) shared(V2,TV2,B,tauV2,U2,TU2,updatetimeB)&
         !$omp& private(starttime,endtime,count_rate,work, ib, jb, nb, istart, istop, istart2, ic)&
         !$omp& firstprivate(ldb,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: B(2,2)) depend(in: V2) depend(in: TV2) depend(in: U2) depend(in: TU2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !
         ! Update B from the right
         !

         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,V2,tauV2,TV2)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )

               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-1)*r)

               if( i+blocksize-1 .ge. istart ) then
                  r1 = max( istart, i )
                  ! Perform the right updates up to the edge of the block
                  do ic = 2,nq
                     istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - nq*r - (nq-ic)*r)
                     r2 = min( i+blocksize-1, istart2-1 )

                     if( r2 .ge. r1) then
                        call dlarf( 'R', r2-r1+1, r, V2( ic, ic, ib ), 1,&
                           tauV2(ic, ib), B( r1, jb+ic ), ldb, work )
                     end if
                  end do
               end if

               r1 = i
               r2 = min( i+blocksize-1, istart-1 )
               call dlarfb( 'R', 'N', 'F', 'C', r2-r1+1, nb, nq, V2(1,1,ib), r+nq, TV2(1,1,ib),nq, &
                  B(r1, jb+1), ldb, work, r2-r1+1 )
            end do
         end do

         !
         ! Update B from the left
         !
         !$omp taskloop default(none) private(istop,c1,c2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,ldb,nblocks,blocksize) &
         !$omp& shared(B,U2,TU2)
         do i = 1,n,blocksize
            c2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )

               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )

               if( c2 .ge. c1 ) then
                  call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, min(nq,nb), U2(1,1,ib), r+nq, TU2(1,1,ib),nq, &
                     B(jb+1, c1), ldb, work, c2-c1+1 )
               end if
            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeB = updatetimeB + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task


         !
         ! Update A
         !

         !$omp task default(none) shared(V2,TV2,A,tauV2,U2,TU2,updatetimeA)&
         !$omp& private(starttime,endtime, count_rate, work, ib, jb, nb, istart, istop, istart2, ic)&
         !$omp& firstprivate(lda,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: A(2,2)) depend(in: V2) depend(in: TV2) depend(in: U2) depend(in: TU2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)
         !
         ! Update A from the right
         !

         !$omp taskloop default(none) private(istart,istart2,r1,r2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,V2,tauV2,TV2)
         do i = 1,n,blocksize
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )

               istart = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-1)*r)

               if( i+blocksize-1 .ge. istart ) then
                  r1 = max( istart, i )
                  ! Perform the right updates up to the edge of the block
                  do ic = 2,nq
                     istart2 = max(jblock + 1,jblock + 1 + (ib-1)*r - (2*nq-ic)*r)
                     r2 = min( i+blocksize-1, istart2-1 )

                     if( r2 .ge. r1 ) then
                        call dlarf( 'R', r2-r1+1, r, V2( ic, ic, ib ), 1,&
                           tauV2(ic, ib), A( r1, jb+ic ), lda, work )
                     end if
                  end do
               end if

               r1 = i
               r2 = min( i+blocksize-1, istart-1 )
               call dlarfb( 'R', 'N', 'F', 'C', r2-r1+1, nb, min(nq,nb), V2(1,1,ib), r+nq, TV2(1,1,ib),nq, &
                  A(r1, jb+1), lda, work, r2-r1+1 )
            end do
         end do

         !
         ! Update A from the left
         !
         !$omp taskloop default(none) private(istop,c1,c2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,lda,nblocks,blocksize) &
         !$omp& shared(A,U2,TU2)
         do i = 1,n,blocksize
            c2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1

               istop = min( n, jblock + (ib-1)*r  + r + nq + r*nq - 1 )
               c1 = max( istop+1, i )

               if( c2 .ge. c1 ) then

                  if( ib .eq. 1 ) then
                     jb = jblock
                     nb = min( n-jb, r+nq )
                     call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, nq, U2(1,1,ib), r+nq, TU2(1,1,ib),nq, &
                        A(jb+1, c1), lda, work, c2-c1+1 )
                  else
                     jb = jblock + 1 + (ib-2)*r
                     nb = min( n-jb-r+1, r+nq )
                     call dlarfb( 'L', 'T', 'F', 'C', nb, c2-c1+1, nq, U2(1,1,ib), r+nq, TU2(1,1,ib),nq, &
                        A(jb+r, c1), lda, work, c2-c1+1  )
                  end if
               end if

            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeA = updatetimeA + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task


         !
         ! Update Q
         !
         !$omp task default(none) shared(U2,TU2,Q,updatetimeQ)&
         !$omp& private(i, work, ib, jb, nb, starttime, endtime, count_rate)&
         !$omp& firstprivate(ldq,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: Q) depend(in: U2) depend(in: TU2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)


         !$omp taskloop default(none) private(r1,r2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,ldq,nblocks,blocksize) &
         !$omp& shared(Q,U2,TU2)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq-1 )
               call dlarfb( 'R', 'N', 'F', 'C', r2-r1+1, nb, nq, U2(1,1,ib), r+nq, TU2(1,1,ib),nq, &
                  Q(r1, jb+1), ldq, work, r2-r1+1 )

            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeQ = updatetimeQ + ((endtime-starttime)/real(count_rate,8))
         end if

         !$omp end task

         !
         ! Update Z
         !
         !$omp task default(none) shared(V2,TV2,Z,updatetimeZ)&
         !$omp& private(i, work, ib, jb, nb, starttime, endtime, count_rate)&
         !$omp& firstprivate(ldz,n,r,nq,jblock,nblocks,blocksize)&
         !$omp& depend(inout: Z) depend(in: V2) depend(in: TV2)&
         !$omp& priority(1)

         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none) private(r1,r2,i,ib,jb,ic,work,nb)&
         !$omp& firstprivate(n,jblock,r,nq,ldz,nblocks,blocksize) &
         !$omp& shared(Z,V2,TV2)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( n, i + blocksize-1 )
            do ib = nblocks,1,-1
               jb = jblock + (ib-1)*r
               nb = min( n-jb, r+nq )
               call dlarfb( 'R', 'N', 'F', 'C', r2-r1+1, nb, nq, V2(1,1,ib), r+nq,&
                  TV2(1,1,ib),nq, Z(r1, jb+1), ldz, work, r2-r1+1 )
            end do
         end do

         if(time_tasks) then
            call system_clock(endtime,count_rate)
            updatetimeZ = updatetimeZ + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            applicationtime = applicationtime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

      end do

      !$omp taskwait
      !$omp end single nowait
      !$omp end parallel

      if(time_tasks) then
         write(*,'("gen: ",F16.6," ;lookahead A: ", F16.6," ;lookahead B: ", F16.6,&
         &" ;update A: ", F16.6," ;update B: ", F16.6," ;update Q: ", F16.6," ;update Z: ", F16.6)')&
            generationtime, lookaheadtimeA, lookaheadtimeB, updatetimeA, updatetimeB, updatetimeQ, updatetimeZ
      end if
      if(time_components) then
         write(*,'("      gen: ",F16.6," ;lookahead: ", F16.6," ;application: ", F16.6)')&
            generationtime2, lookaheadtime, applicationtime
      end if

      deallocate( TV, TU, U, V, tauU, tauV, work, B_copy, TV2, TU2, U2, V2 )

   end subroutine

   subroutine rht2ht_householder_unblocked( n, A, lda, B, ldb, Q, ldq, Z, ldz, r )
      implicit none
      ! arguments
      integer, intent(in) :: n, lda, ldb, ldq, ldz, r
      double precision, intent(inout) :: A(lda,*), B(ldb,*), Z(ldz,*),Q(ldq,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0

      ! local variables
      integer :: i, jc, jb, ierr, nb
      double precision :: tau, temp
      double precision, allocatable :: work(:),B_copy(:,:)

      ! ---------------
      ! executable code
      ! ---------------

      allocate(work( 2*r*n ))
      allocate(B_copy( r, r ))

      ! jc indicates the column we are reducing this iteration
      do jc = 1, n-2

         nb = min( r, n-jc )

         !
         ! Reduce a column in A
         !

         ! Reflector to reduce a column in A
         call dlarfg( nb, A( jc+1, jc ), A( jc+2, jc ), 1, tau )
         temp = A( jc+1, jc )
         A( jc+1, jc ) = one
         ! Apply reflector to A
         call dlarf( 'L', nb, n-jc, A( jc+1, jc ), 1,&
            tau, A( jc+1, jc+1 ), lda, work )
         ! Apply reflector to B
         call dlarf( 'L', nb, n-jc, A( jc+1, jc ), 1,&
            tau, B( jc+1, jc+1 ), ldb, work )
         ! Apply reflector to Q
         call dlarf( 'R', n, nb, A( jc+1, jc ), 1,&
            tau, Q( 1, jc+1 ), ldq, work )
         A( jc+1, jc ) = temp
         do i = 2,nb
            A( jc+i, jc ) = zero
         end do

         ! RQ factorization of block in B
         call dlacpy( 'A', nb, nb, B( jc+1, jc+1 ), ldb, B_copy, r )
         call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
         ! Form orthogonal matrix
         call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
         ! Reflector to annihilate top row of orthogonal matrix
         call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tau )
         B_copy( 1, 1 ) = one
         ! Apply reflector to A
         call dlarf( 'R', n, nb, B_copy(1,1), r,&
            tau, A( 1, jc+1 ), lda, work )
         ! Apply reflector to B
         call dlarf( 'R', jc+nb, nb, B_copy(1,1), r,&
            tau, B( 1, jc+1 ), ldb, work )
         ! Apply reflector to Z
         call dlarf( 'R', n, nb, B_copy(1,1), r,&
            tau, Z( 1, jc+1 ), ldz, work )

         do i = 2,nb
            B( jc+i, jc+1 ) = zero
         end do

         !
         ! Chase down the fill in
         !

         ! jb points to the column in A whose fill in is to be reduced in this iteration
         do jb = jc+1, n-r, r

            nb = min( r, n-jb-r+1 )

            ! Reflector to reduce a column in A
            call dlarfg( nb, A( jb+r, jb ), A( jb+r+1, jb ), 1, tau )
            temp = A( jb+r, jb )
            A( jb+r, jb ) = one
            ! Apply reflector to A
            call dlarf( 'L', nb, n-jb, A( jb+r, jb ), 1,&
               tau, A( jb+r, jb+1 ), lda, work )
            ! Apply reflector to B
            call dlarf( 'L', nb, n-jb-r+1, A( jb+r, jb ), 1,&
               tau, B( jb+r, jb+r ), ldb, work )
            ! Apply reflector to Q
            call dlarf( 'R', n, nb, A( jb+r, jb ), 1,&
               tau, Q( 1, jb+r ), ldq, work )
            A( jb+r, jb ) = temp
            do i = 1,nb-1
               A( jb+r+i, jb ) = zero
            end do

            ! RQ factorization of block in B
            call dlacpy( 'A', nb, nb, B( jb+r, jb+r ), ldb, B_copy, r )
            call dgerq2( nb, nb, B_copy, r, work, work( nb + 1 ), ierr )
            ! Form orthogonal matrix
            call dorgr2( nb, nb, nb, B_copy, r, work, work(nb+1), ierr )
            ! Reflector to annihilate top row of orthogonal matrix
            call dlarfg( nb, B_copy(1,1), B_copy( 1,2 ), r, tau )
            B_copy( 1, 1 ) = one
            ! Apply reflector to A
            call dlarf( 'R', n, nb, B_copy(1,1), r,&
               tau, A( 1, jb+r ), lda, work )
            ! Apply reflector to B
            call dlarf( 'R', jb+r+nb-1, nb, B_copy(1,1), r,&
               tau, B( 1, jb+r ), ldb, work )
            ! Apply reflector to Z
            call dlarf( 'R', n, nb, B_copy(1,1), r,&
               tau, Z( 1, jb+r ), ldz, work )

            do i = 1,nb-1
               B( jb+r+i, jb+r ) = zero
            end do

         end do

      end do


      deallocate( work, B_copy )

   end subroutine

end module

