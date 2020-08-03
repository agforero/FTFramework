
 implicit none

 integer :: m,n
 real(8), allocatable :: matr(:,:)

 interface
   subroutine checkmatr(matr)
   real(8) :: matr(:,:)
   end subroutine checkmatr
 end interface

 write(*,'(a)',advance='no')'Give matrix dimensions m,n: '; read*,m,n
 allocate(matr(m,n))
 call checkmatr(matr)

 end

 subroutine checkmatr(matr)

 implicit none

 real(8) :: matr(:,:)
 real(8) :: localmatr(size(matr,1),size(matr,2))

 print*,size(localmatr)
 print*,shape(localmatr)
 print*,size(localmatr,1),size(localmatr,2)

 end subroutine checkmatr


