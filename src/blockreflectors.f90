module blockreflectors

contains

   !
   !  Apply a block reflector using tasks. Split to account for overlap
   !
   subroutine blockreflectorright( m, n, k, k2, V, ldv, T, ldt, C, ldc, Y, ldy, bl )
      implicit none
      ! arguments
      integer, intent(in) :: m,n,k,k2,ldv,ldt,ldc,ldy, bl
      double precision, intent(inout) :: V(ldv,*), T(ldt,*), C(ldc,*),Y(ldy,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0

      !
      ! V is split into V1 (1:k), V2( k+1:k2 ), V3( k2+1:n )
      !

      !
      ! Form Y1 = C1V1 + C2V2 (does not depend on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(out: Y)&
      !$omp& depend(in: C(1:m,1:k2))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y1 start",bl
      call dlacpy( 'A', m, k, C, ldc, Y, ldy)
      call dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,&
         k, one, V, ldv, Y, ldy )
      if( k2 > k ) then
         call dgemm( 'No transpose', 'No transpose', m, k, k2-k,&
            one, C(1, k+1), ldc, V(k+1,1), ldv, one, Y, ldy)
      end if
      ! write(*,*) "Y1 end",bl
      !$omp end task

      !
      ! Form W = TY = T(Y1 + C3V3) (depends on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y,T) firstprivate(ldc,ldv,ldy,ldt,m,n,k,k2,bl)&
      !$omp& depend(inout: Y)&
      !$omp& depend(in: C(1:m,k2+1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y3 start",bl
      if( n > k2 ) then
         call dgemm( 'No transpose', 'No transpose', m, k, n-k2,&
            one, C(1, k2+1), ldc, V(k2+1,1), ldv, one,&
            Y, ldy)
      end if
      call dtrmm( 'Right', 'Upper', 'No transpose', 'Non-unit', m, k,&
         one, T, ldt, Y, ldy)
      ! write(*,*) "Y3 end",bl
      !$omp end task

      !
      ! Form C1 = C1 - W1 V1**T (The next block reflectors depend on this so
      ! we calculate it first)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: Y)&
      !$omp& depend(inout: C(1:m,1:k))&
      !$omp& priority(bl+1)
      ! write(*,*) "C1 start",bl
      ! Use gemm insteqd of dtrmm, this increases the flops, but
      ! it is out of place so we can do it separately from the other multiplication
      call dgemm( 'No transpose', 'Transpose', m, k, k,&
         -one, Y, ldy, V, ldv, one, C, ldc )
      ! write(*,*) "C1 end",bl
      !$omp end task

      !
      ! Form C23 = C23 - W23 V23**T
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: Y)&
      !$omp& depend(inout: C(1:m,k+1:n))&
      !$omp& priority(bl)
      ! write(*,*) "C23 start",bl
      call dgemm( 'No transpose', 'Transpose', m, n-k, k,&
         -one, Y, ldy, V(k+1,1), ldv, one, C(1,k+1), ldc )
      ! write(*,*) "C23 end",bl
      !$omp end task

   end subroutine
   !
   !  Apply a block reflector using tasks. Split to account for overlap
   !
   subroutine blockreflectorrightrb( m, n, k, k2, k3, V, ldv, T, ldt, C, ldc, Y, ldy, bl )
      implicit none
      ! arguments
      integer, intent(in) :: m,n,k,k2,k3,ldv,ldt,ldc,ldy, bl
      double precision, intent(inout) :: V(ldv,*), T(ldt,*), C(ldc,*),Y(ldy,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0

      integer :: i, j

      !
      ! V is split into V1 (1:k), V2( k+1:k2 ), V3( k2+1:n )
      !

      !
      ! Form Y1 = C1 V1**T + C2 V2**T (does not depend on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,k3,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(out: Y)&
      !$omp& depend(in: C(1:k,1:k2))&
      !$omp& depend(in: C(k+1:k3,1:k2))&
      !$omp& depend(in: C(k3+1:m,1:k2))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y1 start",bl
      call dgemm( 'No transpose', 'Transpose', m, k, k2,&
         one, C, ldc, V, ldv, zero, Y, ldy)
      ! write(*,*) "Y1 end",bl
      !$omp end task

      !
      ! Form W = YT = (Y1 + C3 V3**T)T (depends on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y,T) firstprivate(ldc,ldv,ldy,ldt,m,n,k,k2,k3,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(inout: Y)&
      !$omp& depend(in: C(1:k,k2+1:n))&
      !$omp& depend(in: C(k+1:k3,k2+1:n))&
      !$omp& depend(in: C(k3+1:m,k2+1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y3 start",bl
      ! if( n > k2 ) then
      !    call dgemm( 'No transpose', 'No transpose', m, k, n-k2,&
      !       one, C(1, k2+1), ldc, V(k2+1,1), ldv, one,&
      !       Y, ldy)
      ! end if
      if( n > k2 ) then
         call dgemm( 'No transpose', 'Transpose', m, k, n-k2,&
            one, C(1,k2+1), ldc, V(1,k2+1), ldv, one, Y, ldy)
      end if
      call dtrmm( 'Right', 'Lower', 'No transpose', 'Non-unit', m, k,&
         one, T, ldt, Y, ldy)
      ! write(*,*) "Y3 end",bl
      !$omp end task

      !
      ! Form C = C - W V
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,k3,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(inout: Y)&
      !$omp& depend(inout: C(1:k,1:k))&
      !$omp& depend(inout: C(1:k,k+1:n))&
      !$omp& depend(inout: C(k+1:k3,1:k))&
      !$omp& depend(inout: C(k+1:k3,k+1:n))&
      !$omp& depend(inout: C(k3+1:m,1:k))&
      !$omp& depend(inout: C(k3+1:m,k+1:n))&
      !$omp& priority(bl)
      call dgemm( 'No transpose', 'No transpose', m, n-k, k,&
         -one, Y, ldy, V, ldv, one, C, ldc )
      call dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,&
         k, one, V( 1, n-k+1 ), ldv, Y, ldy )
      do j = 1,k
         do i=1,m
            C(i,n-k+j) = C(i,n-k+j) - Y(i,j)
         end do
      end do
      !$omp end task

   end subroutine

   !
   !  Apply a block reflector using tasks. Split to account for overlap
   !
   subroutine blockreflectorrightrf( m, n, k, k2, k3, V, ldv, T, ldt, C, ldc, Y, ldy, bl )
      implicit none
      ! arguments
      integer, intent(in) :: m,n,k,k2,ldv,ldt,ldc,ldy, bl, k3
      double precision, intent(inout) :: V(ldv,*), T(ldt,*), C(ldc,*),Y(ldy,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0

      integer :: i,j

      !
      ! V is split into V1 (1:k), V2( k+1:k2 ), V3( k2+1:n )
      !

      !
      ! Form Y1 = C1 V1**T + C2 V2**T (does not depend on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(out: Y)&
      !$omp& depend(in: C(1:k3,1:k2))&
      !$omp& depend(in: C(k3+1:m,1:k2))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y1 start",bl
      do j=1,k
         do i=1,m
            Y(i,j) = C(i,j)
         end do
      end do
      call dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', m,&
         k, one, V, ldv, Y, ldy )
      if( k2 > k ) then
         call dgemm( 'No transpose', 'Transpose', m, k, k2-k,&
            one, C(1, k+1), ldc, V(1,k+1), ldv, one, Y, ldy)
      end if
      ! write(*,*) "Y1 end",bl
      !$omp end task

      !
      ! Form W = TY = (Y1 + C3 V3**T)T (depends on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y,T) firstprivate(ldc,ldv,ldy,ldt,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(inout: Y)&
      !$omp& depend(in: C(1:k3,k2+1:n))&
      !$omp& depend(in: C(k3+1:m,k2+1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y3 start",bl
      if( n > k2 ) then
         call dgemm( 'No transpose', 'Transpose', m, k, n-k2,&
            one, C(1, k2+1), ldc, V(1,k2+1), ldv, one,&
            Y, ldy)
      end if
      call dtrmm( 'Right', 'Upper', 'No transpose', 'Non-unit', m, k,&
         one, T, ldt, Y, ldy)
      ! write(*,*) "Y3 end",bl
      !$omp end task

      !
      ! Form C1 = C1 - W1 V1 (The next block reflectors depend on this so
      ! we calculate it first)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(in: Y)&
      !$omp& depend(inout: C(1:k3,1:k))&
      !$omp& depend(inout: C(k3+1:m,1:k))&
      !$omp& priority(bl+1)
      ! write(*,*) "C1 start",bl
      ! Use gemm insteqd of dtrmm, this increases the flops, but
      ! it is out of place so we can do it separately from the other multiplication
      call dgemm( 'No transpose', 'No transpose', m, k, k,&
         -one, Y, ldy, V, ldv, one, C, ldc )
      ! write(*,*) "C1 end",bl
      !$omp end task

      !
      ! Form C23 = C23 - W23 V23**T
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(in: Y)&
      !$omp& depend(inout: C(1:k3,k+1:n))&
      !$omp& depend(inout: C(k3+1:m,k+1:n))&
      !$omp& priority(bl)
      ! write(*,*) "C23 start",bl
      call dgemm( 'No transpose', 'No transpose', m, n-k, k,&
         -one, Y, ldy, V(1,k+1), ldv, one, C(1,k+1), ldc )
      ! write(*,*) "C23 end",bl
      !$omp end task

   end subroutine



   !
   !  Apply a block reflector using tasks. Split to account for overlap
   !
   subroutine blockreflectorleft( m, n, k, k2, V, ldv, T, ldt, C, ldc, Y, ldy, bl )
      implicit none
      ! arguments
      integer, intent(in) :: m,n,k,k2,ldv,ldt,ldc,ldy, bl
      double precision, intent(inout) :: V(ldv,*), T(ldt,*), C(ldc,*),Y(ldy,*)

      ! parameters
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0

      integer :: i,j

      !
      ! V is split into V1 (1:k), V2( k+1:k2 ), V3( k2+1:n )
      !

      !
      ! Form Y1 = C1**TV1 + C2**TV2 (does not depend on previous block reflectors)
      !

      !$omp task default(none) private(i,j) shared(C,V,Y)&
      !$omp& firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(out: Y)&
      !$omp& depend(in: C(1:k2,1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y1 start",bl
      do j=1,k
         do i=1,n
            Y(i,j) = C(j,i)
         end do
      end do
      call dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,&
         k, one, V, ldv, Y, ldy )
      if( k2 > k ) then
         call dgemm( 'Transpose', 'No transpose', n, k, k2-k,&
            one, C(k+1, 1), ldc, V(k+1,1), ldv, one, Y, ldy)
      end if
      ! write(*,*) "Y1 end",bl
      !$omp end task

      !
      ! Form W = YT**T = (Y1 + C3**TV3)T**T (depends on previous block reflectors)
      !

      !$omp task default(none) shared(C,V,Y,T) firstprivate(ldc,ldv,ldy,ldt,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(inout: Y)&
      !$omp& depend(in: C(k2+1:m,1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "Y3 start",bl
      if( m > k2 ) then
         call dgemm( 'Transpose', 'No transpose', n, k, m-k2,&
            one, C(k2+1,1), ldc, V(k2+1,1), ldv, one,&
            Y, ldy)
      end if
      call dtrmm( 'Right', 'Upper', 'No transpose', 'Non-unit', n, k,&
         one, T, ldt, Y, ldy)
      ! write(*,*) "Y3 end",bl
      !$omp end task

      !
      ! Form C1 = C1 - V1 W1**T (The next block reflectors depend on this so
      ! we calculate it first)
      !

      !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
      !$omp& depend(in: V, T)&
      !$omp& depend(in: Y)&
      !$omp& depend(inout: C(1:k,1:n))&
      !$omp& priority(bl+1)
      ! write(*,*) "C1 start",bl
      ! Use gemm insteqd of dtrmm, this increases the flops, but
      ! it is out of place so we can do it separately from the other multiplication
      call dgemm( 'No transpose', 'Transpose', k, n, k,&
         -one, V, ldv, Y, ldy, one, C, ldc )
      ! write(*,*) "C1 end",bl
      !$omp end task

      !
      ! Form C23 = C23 - V23 W23**T
      !

      if(m > k) then
         !$omp task default(none) shared(C,V,Y) firstprivate(ldc,ldv,ldy,m,n,k,k2,bl)&
         !$omp& depend(in: V, T)&
         !$omp& depend(in: Y)&
         !$omp& depend(inout: C(k+1:k2,1:n))&
         !$omp& depend(inout: C(k2+1:m,1:n))&
         !$omp& priority(bl)
         ! write(*,*) "C23 start",bl
         call dgemm( 'No transpose', 'Transpose', m-k, n, k,&
            -one, V(k+1,1), ldv, Y, ldy, one, C(k+1,1), ldc )
         ! write(*,*) "C23 end",bl
         !$omp end task
      end if

   end subroutine



end module
