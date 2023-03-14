module ab2rht
   use omp_lib
   use utils
   use blockreflectors

contains

   subroutine ab2rht_simple( n, A, lda, B, ldb, Q, ldq, Z, ldz, nb, p, time_components )
      implicit none
      ! arguments
      integer, intent(in) :: n, lda, ldb, ldq, ldz, nb, p
      double precision, intent(inout) :: A(lda,*), B(ldb,*), Z(ldz,*),Q(ldq,*)
      logical :: time_components

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, parameter :: min_blocksize = 32

      ! local variables
      integer :: i,j, ib, jb, nblocks, ib2, ierr, nb2, lwork, k, bl,r1,r2,c1,c2, nthreads, blocksize
      double precision, allocatable :: VL(:,:,:), VR(:,:,:), TL(:,:,:), TR(:,:,:), work(:), B_copy(:,:)
      double precision :: dummywork(1)

      ! temporary timing variables
      logical, parameter :: time_tasks = .false.
      integer(8) :: starttime,endtime,count_rate
      double precision :: genlefttime, genrighttime, leftmulAtime, leftmulBtime,&
         mulQtime, rightmulAtime, mulZtime, multime, genlefttime2, genrighttime2

      ! ---------------
      ! executable code
      ! ---------------
      nblocks = (n-nb)/((p-1)*nb)
      if( mod( (n-nb) , ((p-1)*nb) ) .ne. 0 ) nblocks = nblocks + 1
      allocate( VL( p*nb,nb, nblocks ) )
      allocate( VR( nb,p*nb, nblocks ) )
      allocate( TL( nb, nb, nblocks ) )
      allocate( TR( nb, nb, nblocks ) )
      allocate( B_copy( p*nb, p*nb ) )

      lwork = p*nb*n
      call dgerqf( p*nb, p*nb, B_copy, p*nb, dummywork, dummywork, -1, ierr )
      lwork = max( lwork, int(dummywork(1)) )
      call dorgrq( p*nb, p*nb, p*nb, B_copy, p*nb, dummywork, dummywork, -1, ierr )
      lwork = max( lwork, int(dummywork(1)) )
      allocate( work( lwork ) )

      genlefttime = 0.0d0
      genrighttime = 0.0d0
      leftmulAtime = 0.0d0
      leftmulBtime = 0.0d0
      mulQtime = 0.0d0
      rightmulAtime = 0.0d0
      mulZtime = 0.0d0
      multime = 0.0d0
      genlefttime2 = 0.0d0
      genrighttime2 = 0.0d0

      !$omp parallel
      !$omp single

      nthreads = omp_get_max_threads()
      blocksize = n/nthreads
      if( mod( blocksize, 16 ) .ne. 0 ) blocksize = blocksize + (16-mod( blocksize, 16 ))
      blocksize = max( min_blocksize, blocksize )
      blocksize = min( n, blocksize )


      ! do jb = 1, 1, nb
      do jb = 1, n-nb-1, nb
         ! Reduce panel A(jb+nb:n,jb:jb+nb-1)

         nblocks = (n-nb-jb+1)/((p-1)*nb)
         if( mod( (n-nb-jb+1) , ((p-1)*nb) ) .ne. 0 ) nblocks = nblocks + 1

         if( time_components ) then
            !$omp taskwait
            call system_clock(starttime)
         end if

         !
         ! Generate multiplications from the left
         !
         !$omp task default(none)&
         !$omp& shared(A,TL,VL,genlefttime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,lda)&
         !$omp& private(bl,ib,ib2,nb2,k,work,ierr,starttime,endtime,count_rate)&
         !$omp& depend(inout: A)&
         !$omp& depend(out: VL) depend(out: TL)&
         !$omp& priority(100)
         if(time_tasks) call system_clock(starttime)
         bl = 0
         do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

            bl = bl + 1
            ib2 = min( n, ib + p*nb - 1 )
            nb2 = ib2-ib+1
            k = min(nb, nb2)

            !
            ! Generate left reflectors to reduce block of panel
            ! QR of A( ib:ib2, jb:jb+nb-1 )
            !
            call dgeqr2( nb2, nb, A( ib, jb ), lda, work, work(nb+1), ierr )

            ! Copy reflectors in A to V so we can use them in parallel
            call dlacpy( 'A', nb2, nb, A(ib,jb),lda, VL(1,1,bl), p*nb )
            do j = jb,jb+nb-1
               do i = ib + j-jb+1,ib2
                  A(i,j) = 0.0d0
               end do
            end do

            ! Form T
            call dlarft( 'F', 'C', nb2, k, VL(1,1,bl), p*nb, work, TL(1,1,bl), nb )

         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            genlefttime = genlefttime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            genlefttime2 = genlefttime2 + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

         ! 
         ! Apply left multiplications to B
         !
         !$omp task default(none)&
         !$omp& shared(B,TL,VL,leftmulBtime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldb,blocksize)&
         !$omp& private(c1,c2,i,bl,ib,ib2,nb2,k,work,starttime,endtime,count_rate)&
         !$omp& depend(inout: B)&
         !$omp& depend(in: VL) depend(in: TL)&
         !$omp& priority(100)
         if(time_tasks) call system_clock(starttime)


         !$omp taskloop default(none)&
         !$omp& shared(B,TL,VL)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldb,blocksize)&
         !$omp& private(c1,c2,i,bl,ib,ib2,nb2,k,work)
         do i = 1,n,blocksize
            c2 = min( i+blocksize-1, n )
            bl = 0
            do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

               bl = bl + 1
               ib2 = min( n, ib + p*nb - 1 )
               nb2 = ib2-ib+1
               k = min(nb, nb2)
               c1 = max(ib,i)

               if( c2 .ge. c1 ) then
                  call dlarfb( 'L', 'T', 'F', 'C', nb2, c2-c1+1, k, &
                     VL(1,1,bl), p*nb, TL(1,1,bl), nb, B(ib, c1), ldb, work, c2-c1+1 )
               end if

            end do
         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            leftmulBtime = leftmulBtime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         ! 
         ! Apply left multiplications to A
         ! 
         !$omp task default(none)&
         !$omp& shared(A,TL,VL,leftmulAtime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,lda,blocksize)&
         !$omp& private(bl,ib,ib2,nb2,k,work,starttime,endtime,count_rate)&
         !$omp& depend(inout: A)&
         !$omp& depend(in: VL) depend(in: TL)&
         !$omp& priority(1)
         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none)&
         !$omp& shared(A,TL,VL)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,lda,blocksize)&
         !$omp& private(c1,c2,i,bl,ib,ib2,nb2,k,work)
         do i = 1,n,blocksize
            c1 = max( i, jb+nb )
            c2 = min( i+blocksize-1, n )

            if( c2 .ge. c1 ) then
               bl = 0
               do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

                  bl = bl + 1
                  ib2 = min( n, ib + p*nb - 1 )
                  nb2 = ib2-ib+1
                  k = min(nb, nb2)

                  call dlarfb( 'L', 'T', 'F', 'C', nb2, c2-c1+1, k, &
                     VL(1,1,bl), p*nb, TL(1,1,bl), nb, A(ib, c1), lda, work, c2-c1+1 )

               end do
            end if
         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            leftmulAtime = leftmulAtime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         ! 
         ! Apply left multiplications to Q
         ! 
         !$omp task default(none)&
         !$omp& shared(Q,TL,VL,mulQtime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldq,blocksize)&
         !$omp& private(bl,ib,ib2,nb2,k,work,starttime,endtime,count_rate)&
         !$omp& depend(inout: Q)&
         !$omp& depend(in: VL) depend(in: TL)&
         !$omp& priority(1)
         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none)&
         !$omp& shared(Q,TL,VL)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldq, blocksize)&
         !$omp& private(r1,r2,i,bl,ib,ib2,nb2,k,work)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( i+blocksize-1, n )
            bl = 0
            do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

               bl = bl + 1
               ib2 = min( n, ib + p*nb - 1 )
               nb2 = ib2-ib+1
               k = min(nb, nb2)

               call dlarfb( 'R', 'N', 'F', 'C', r2-r1+1, nb2, k, &
                  VL(1,1,bl), p*nb, TL(1,1,bl), nb, Q(r1, ib), ldq, work, r2-r1+1 )

            end do
         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            mulQtime = mulQtime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            multime = multime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

         !
         ! Generate multiplications from the right to reduce fill in
         ! In the process, also applies these updates to B
         !
         !$omp task default(none)&
         !$omp& shared(B,TR,VR,B_copy,genrighttime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldb,lwork,blocksize)&
         !$omp& private(bl,ib,ib2,nb2,k,work,ierr,starttime,endtime,count_rate)&
         !$omp& depend(inout: B)&
         !$omp& depend(out: VR) depend(out: TR)&
         !$omp& priority(100)
         if(time_tasks) call system_clock(starttime)
         bl = 0
         do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

            bl = bl + 1
            ib2 = min( n, ib + p*nb - 1 )
            nb2 = ib2-ib+1
            k = min(nb, nb2)

            ! 
            ! Actual generation of the reflectors
            ! 
            !$omp task default(none)&
            !$omp& shared(B,TR,VR,B_copy)&
            !$omp& firstprivate(bl,ib,ib2,nb2,k,jb,n,nb,p,nblocks,ldb,lwork)&
            !$omp& private(work,ierr)&
            !$omp& depend(in: B(2,2))&
            !$omp& depend(out: VR(1,1,bl)) depend(out: TR(1,1,bl))

            ! RQ factorization of subblock of B to remove fill in
            call dlacpy( 'A', nb2, nb2, B(ib,ib), ldb, B_copy, p*nb )
            call dgerqf( nb2, nb2, B_copy, p*nb, work, work(nb2+1), lwork-nb2, ierr )

            if( nb2 .gt. nb ) then

               ! LQ factorization of orthogonal factor in preceding RQ
               ! only reduces nb instead of nb2 columns (much cheaper)
               call dorgrq( nb2, nb2, nb2, B_copy, p*nb, work, work(p*nb+1), lwork-nb2, ierr )
               call dgelq2( nb, nb2, B_copy, p*nb, work, work(nb+1), ierr )
               call dlacpy( 'A', nb, nb2, B_copy, p*nb, VR(1,1,bl), nb )

               ! Form T
               call dlarft( 'F', 'R', nb2, nb, VR(1,1,bl), nb, work, TR(1,1,bl), nb )

            else
               ! nb2 is small, just apply the RQ factorization without further tricks

               ! Form T
               call dlarft( 'B', 'R', nb2, nb2, B_copy, p*nb, work, TR(1,1,bl), nb )
               call dlacpy( 'A', nb2, nb2, B_copy, p*nb, VR(1,1,bl), nb )

            end if
            !$omp end task

            ! 
            ! Lookahead update of B
            ! 
            !$omp task default(none)&
            !$omp& shared(B,TR,VR)&
            !$omp& firstprivate(bl,ib,ib2,nb2,k,jb,n,nb,p,nblocks,ldb)&
            !$omp& private(c1,c2,i,j,work)&
            !$omp& depend(inout: B(2,2)) depend(in: B(3,3))&
            !$omp& depend(in: VR(1,1,bl)) depend(in: TR(1,1,bl))

            c1 = max( 1, ib-nb*(p-1) )
            c2 = ib2
            if( nb2 .gt. nb ) then
               call dlarfb( 'R', 'N', 'F', 'R', c2-c1+1, nb2, nb, &
                  VR(1,1,bl), nb, TR(1,1,bl), nb, B(c1, ib), ldb, work, c2-c1+1 )
               do j = ib,ib+nb-1
                  do i = j+1,ib2
                     B(i,j) = 0.0d0
                  end do
               end do
            else
               call dlarfb( 'R', 'N', 'B', 'R', c2-c1+1, nb2, nb2, &
                  VR(1,1,bl), nb, TR(1,1,bl), nb, B(c1, ib), ldb, work, c2-c1+1 )
            end if
            !$omp end task

            ! 
            ! Update the rest of B
            ! 
            !$omp task default(none)&
            !$omp& shared(B,TR,VR)&
            !$omp& firstprivate(bl,ib,ib2,nb2,k,jb,n,nb,p,nblocks,ldb,blocksize)&
            !$omp& private(c1,c2,work)&
            !$omp& depend(inout: B(3,3))&
            !$omp& depend(in: VR(1,1,bl)) depend(in: TR(1,1,bl))

            !$omp taskloop default(none)&
            !$omp& shared(B,TR,VR)&
            !$omp& firstprivate(bl,ib,ib2,nb2,k,jb,n,nb,p,nblocks,ldb,blocksize)&
            !$omp& private(c1,c2,i,work)
            do i = 1,ib,blocksize
               c1 = i
               c2 = min( i+blocksize-1, max( 1, ib-nb*(p-1) )-1)
               if(c2 .ge. c1) then
                  if( nb2 .gt. nb ) then
                     call dlarfb( 'R', 'N', 'F', 'R', c2-c1+1, nb2, nb, &
                        VR(1,1,bl), nb, TR(1,1,bl), nb, B(c1, ib), ldb, work, c2-c1+1 )
                  else
                     call dlarfb( 'R', 'N', 'B', 'R', c2-c1+1, nb2, nb2, &
                        VR(1,1,bl), nb, TR(1,1,bl), nb, B(c1, ib), ldb, work, c2-c1+1 )
                  end if
               end if
            end do
            !$omp end task

         end do
         !$omp taskwait
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            genrighttime = genrighttime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            genrighttime2 = genrighttime2 + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

         !
         ! Apply the right updates to A
         ! 
         !$omp task default(none)&
         !$omp& shared(A,TR,VR,rightmulAtime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,lda,blocksize)&
         !$omp& private(bl,ib,ib2,nb2,k,work,starttime,endtime,count_rate)&
         !$omp& depend(inout: A)&
         !$omp& depend(in: VR) depend(in: TR)&
         !$omp& priority(1)
         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none)&
         !$omp& shared(A,TR,VR)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,lda,blocksize)&
         !$omp& private(r1,r2,i,bl,ib,ib2,nb2,k,work)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( i+blocksize-1, n )
            bl = 0
            do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

               bl = bl + 1
               ib2 = min( n, ib + p*nb - 1 )
               nb2 = ib2-ib+1
               k = min(nb, nb2)

               if( nb2 .gt. nb ) then
                  call dlarfb( 'R', 'N', 'F', 'R', r2-r1+1, nb2, nb, &
                     VR(1,1,bl), nb, TR(1,1,bl), nb, A(r1, ib), lda, work, r2-r1+1 )
               else
                  call dlarfb( 'R', 'N', 'B', 'R', r2-r1+1, nb2, nb2, &
                     VR(1,1,bl), nb, TR(1,1,bl), nb, A(r1, ib), lda, work, r2-r1+1 )
               end if

            end do
         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            rightmulAtime = rightmulAtime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         !
         ! Apply the right updates to Z
         ! 
         !$omp task default(none)&
         !$omp& shared(Z,TR,VR,mulZtime)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldz,blocksize)&
         !$omp& private(bl,ib,ib2,nb2,k,work,starttime,endtime,count_rate)&
         !$omp& depend(inout: Z)&
         !$omp& depend(in: VR) depend(in: TR)&
         !$omp& priority(1)
         if(time_tasks) call system_clock(starttime)

         !$omp taskloop default(none)&
         !$omp& shared(Z,TR,VR)&
         !$omp& firstprivate(jb,n,nb,p,nblocks,ldz,blocksize)&
         !$omp& private(r1,r2,i,bl,ib,ib2,nb2,k,work)
         do i = 1,n,blocksize
            r1 = i
            r2 = min( i+blocksize-1, n )
            bl = 0
            do ib = jb + nb + (nblocks-1)*(p-1)*nb, jb+nb, -nb*(p-1)

               bl = bl + 1
               ib2 = min( n, ib + p*nb - 1 )
               nb2 = ib2-ib+1
               k = min(nb, nb2)

               if( nb2 .gt. nb ) then
                  call dlarfb( 'R', 'N', 'F', 'R', r2-r1+1, nb2, nb, &
                     VR(1,1,bl), nb, TR(1,1,bl), nb, Z(r1, ib), ldz, work, r2-r1+1 )
               else
                  call dlarfb( 'R', 'N', 'B', 'R', r2-r1+1, nb2, nb2, &
                     VR(1,1,bl), nb, TR(1,1,bl), nb, Z(r1, ib), ldz, work, r2-r1+1 )
               end if
            end do
         end do
         if(time_tasks) then
            call system_clock(endtime,count_rate)
            mulZtime = mulZtime + ((endtime-starttime)/real(count_rate,8))
         end if
         !$omp end task

         if( time_components ) then
            !$omp taskwait
            call system_clock(endtime,count_rate)
            multime = multime + ((endtime-starttime)/real(count_rate,8))
            call system_clock(starttime)
         end if

      end do

      !$omp taskwait
      !$omp end single nowait
      !$omp end parallel

      if(time_tasks) then
         write(*,'("      genlefttime: ",F16.6," ;genrighttime: ", F16.6," ;leftmulAtime: ", F16.6,&
            &" ;leftmulBtime: ", F16.6," ;rightmulAtime: ", F16.6," ;mulQtime: ", F16.6," ;mulZtime: ", F16.6)')&
            genlefttime, genrighttime, leftmulAtime, leftmulBtime, rightmulAtime, mulQtime, mulZtime
      end if
      if(time_components) then
         write(*,'("      genlefttime: ",F16.6," ;genrighttime: ", F16.6," ;multime: ", F16.6)')&
            genlefttime2, genrighttime2, multime
      end if


      deallocate( VL, VR, TL, TR, B_copy, work )

   end subroutine

end module

