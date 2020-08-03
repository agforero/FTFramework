
 implicit none

 integer :: f,n
 real(8) :: p,x1,x2,s
 real(8), external :: power,exponential

 1 print*,'Function (1=power, 2=exponential):'; read*,f; 
 if (f/=1.and.f/=2) goto 1
 print*,'Give power p (in x**p or exp(p*x)):'; read*,p
 print*,'Integrate between values (x1,x2):'; read*,x1,x2
 print*,'Number n of integration intervals:'; read*,n
 if (mod(n,2)/=0) then
   n=n+1; print*,'Using n = ',n
 end if
 if (f==1) then
    call simpson(x1,x2,n,power,p,s)
 else
    call simpson(x1,x2,n,exponential,p,s)
 endif
 print*,'The Simpson integral is: ',s

 end

!------------------------------------!
 subroutine simpson(x1,x2,n,func,p,s)

 implicit none
 
 integer :: i,n
 real(8) :: x1,x2,dx,p,s
 real(8), external :: func
 
 dx=(x2-x1)/dble(n)
 s=func(x1,p)
 do i=1,n-1,2
   s=s+func(x1+dx*i,p)*4.d0+func(x1+dx*(i+1),p)*2.d0
 end do
 s=(s-func(x2,p))/(3.d0*n)

 end subroutine simpson

!-------------------!
 function power(x,p)
 implicit none
 real(8) :: power,x,p
 power=x**p
 end function power

!-------------------------!
 function exponential(x,p)
 implicit none
 real(8) :: exponential,x,p
 exponential=exp(p*x)
 end function exponential

