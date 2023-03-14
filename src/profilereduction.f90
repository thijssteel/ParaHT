program test
   use utils
   use ab2rht
   use rht2ht
   USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
   implicit none

   integer :: n, r, p, nq, algorithm, nb
   logical :: time_components

   double precision, allocatable :: A(:,:),B(:,:),Q(:,:),Z(:,:),A_copy(:,:),B_copy(:,:), C(:,:)

   double precision, external :: dlange
   double precision :: anorm, aerror, bnorm, berror, rhttime, httime, totaltime

   integer :: i, j, info
   integer(8) :: starttime,endtime,count_rate

   character(len = 32) :: arg

   integer, allocatable :: seed(:)
   integer :: seed_size

   ! Parse command line arguments
   if( command_argument_count() .lt. 2 ) then
      write(*,*) 'usage: command n algorithm (r p nq time_components)'
      stop
   end if

   call getarg( 1, arg )
   read( arg, * ) n

   call getarg( 2, arg )
   read( arg, * ) algorithm
   ! Algorithm = 1 for ParaHT
   ! Algorithm = 2 for LAPACK

   call random_seed( size = seed_size )
   allocate( seed( seed_size ) )
   do i = 1, seed_size
      seed( i ) = 0
   end do
   seed( 1 ) = 1302
   call random_seed( put = seed )

   if( algorithm .eq. 1 ) then

      if( command_argument_count() .lt. 5 ) then
         write(*,*) 'if using ParaHT, you should supply the parameters r, p and nq'
         stop
      end if

      call getarg( 3, arg )
      read( arg, * ) r

      call getarg( 4, arg )
      read( arg, * ) p

      call getarg( 5, arg )
      read( arg, * ) nq

      time_components = .false.
      if( command_argument_count() .ge. 6 ) then
         call getarg( 6, arg )
         read( arg, * ) i
         if( i .eq. 1 ) time_components = .true.
      end if

   end if

   ! Prepare the matrices

   allocate(A(n,n))
   allocate(B(n,n))
   allocate(Q(n,n))
   allocate(Z(n,n))
   allocate(A_copy(n,n))
   allocate(B_copy(n,n))
   allocate(C(n,n))

   !
   ! Initialize the data in parallel to avoid all the data being
   ! allocated close to one node
   !

   !$omp parallel do
   do j=1,n
      do i=1,n
         A(i,j) = 0.0d0
      end do
   end do

   !$omp parallel do
   do j=1,n
      do i=1,n
         B(i,j) = 0.0d0
      end do
   end do

   !$omp parallel do
   do j=1,n
      do i=1,n
         Q(i,j) = 0.0d0
      end do
   end do

   !$omp parallel do
   do j=1,n
      do i=1,n
         Z(i,j) = 0.0d0
      end do
   end do

   call random_number( A )
   call random_number( B )

   call dgeqrf( n, n, B, n, A_copy, C, n*n, info )

   do j=1,n
      do i=j+1,n
         B(i,j) = 0.0d0
      end do
   end do

   call dlaset( 'A', n, n, 0.0d0, 1.0d0, Q, n )
   call dlaset( 'A', n, n, 0.0d0, 1.0d0, Z, n )

   call dlacpy('A', n, n, A, n, A_copy, n)
   call dlacpy('A', n, n, B, n, B_copy, n)

   ! Timings
   if(algorithm .eq. 1) then

      call system_clock(starttime)

      call ab2rht_simple( n, A, n, B, n, Q, n, Z, n, r, p, time_components )

      call system_clock(endtime,count_rate)

      rhttime = ((endtime-starttime)/real(count_rate,8))

      call system_clock(starttime)

      call rht2ht_householder( n, A, n, B, n, Q, n, Z, n, r, nq, time_components )

      call system_clock(endtime,count_rate)

      httime = ((endtime-starttime)/real(count_rate,8))

      totaltime = rhttime + httime

   else if( algorithm .eq. 2 ) then

      call system_clock(starttime)

      call dgghd3( 'V', 'V', n, 1, n, A, n, B, n, Q, n, Z, n, C, n*n, info )

      call system_clock(endtime,count_rate)

      totaltime = ((endtime-starttime)/real(count_rate,8))

   end if

   ! Check error
   call dgemm( 'T', 'N', n, n, n, 1.0d0, Q, n, A_copy, n, 0.0d0, C, n )
   call dgemm( 'N', 'N', n, n, n, 1.0d0, C, n, Z, n, -1.0d0, A, n )

   call dgemm( 'T', 'N', n, n, n, 1.0d0, Q, n, B_copy, n, 0.0d0, C, n )
   call dgemm( 'N', 'N', n, n, n, 1.0d0, C, n, Z, n, -1.0d0, B, n )

   anorm = dlange( 'F', n, n, A_copy, n, C )
   bnorm = dlange( 'F', n, n, B_copy, n, C )
   aerror = dlange( 'F', n, n, A, n, C )
   berror = dlange( 'F', n, n, B, n, C )

   if( algorithm .eq. 1 ) then

      write(*,'("n:",I6," ;runtime:",F16.6," ;error:",2ES16.6, " ;phase 1 time:",F16.6, " ;phase 2 time:",F16.6,&
         &" ;r:",I4," ;p:",I4," ;nq:",I4)') &
         n, totaltime, aerror/anorm, berror/bnorm,rhttime, httime, r, p, nq
   else

      write(*,'("n:",I6," ;runtime:",F16.6," ;error:",2ES16.6)') &
         n, totaltime, aerror/anorm, berror/bnorm
   end if

   deallocate( A, B, Q, Z, A_copy, B_copy, C, seed )

contains

end program test
