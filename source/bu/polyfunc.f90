
 implicit none

 integer, parameter :: nmax=10

 integer :: i,n
 real(8) :: a(0:nmax),x

 real(8), external :: poly

 print*,'Order of the polynomial'; read*,n
 if (n > nmax) then
   print*,'Order higher than nmax'; stop
 endif

 do i=0,n
   print*,'Give coefficient ',i; read*,a(i)
 enddo

 print*,'Evaluate at x-value:'; read*,x
 print*,'Polynomial value is: ',poly(n,a(0:n),x)

 end

 function poly(n,a,x)

 implicit none

 integer :: i,n
 real(8) :: poly,a(0:n),x
 
 poly=0.0d0
 do i=0,n
   poly=poly+a(i)*x**i
 enddo

 end function poly
