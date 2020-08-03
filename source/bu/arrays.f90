
 implicit none

 integer, dimension(2,2) :: a,b
 integer                 :: c(2)

 a(1,1)=1;  a(2,1)=2;  a(1,2)=3;  a(2,2)=4
 b=2*a+1
 print*,a
 print*,b

 c=a(1,1:2)
 print*,c

 end
