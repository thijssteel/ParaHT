module utils

contains

   subroutine print_matrix(m,n, A, lda)
      implicit none
      integer, intent(in) :: m, n, lda
      double precision, intent(in) :: A(lda, *)

      integer ii,jj

      do ii = 1,m
         do jj = 1,n
            write(*,'(F16.6)',advance='no') A(ii,jj)
         end do
         write(*,*)
      end do

   end subroutine
   
end module

