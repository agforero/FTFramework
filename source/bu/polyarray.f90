module powermodule

  integer, parameter :: nd=10, np=3
  real(8) :: xpowers(np,nd)

contains

  subroutine initpowers(dx)

  implicit none

  integer :: i,j
  real(8) :: dx

  do i=1,nd
     xpowers(1,i)=dx*i
  enddo
  do i=2,np
     xpowers(i,:)=xpowers(1,:)**i
  enddo
   
  end subroutine initpowers

end module powermodule


 use powermodule

 implicit none

 integer :: i
 real(8) :: a(0:np),px,dx

 print*,'Give x-spacing '; read*,dx
 do i=0,np
   print*,'Give polynomial coefficient ',i; read*,a(i)
 enddo
 call initpowers(dx)
 do i=1,nd
   px=a(0)+dot_product(xpowers(1:np,i),a(1:np))
   print*,xpowers(1,i),px
 enddo

 end
 