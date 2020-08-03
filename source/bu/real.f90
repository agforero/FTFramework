 implicit none
  
 real(4) :: a
 real(8) :: b
 real(16) :: c

 print*,kind(1.),kind(1.d0)
 a=3.14159265358979328e0
 print*,a
 b=3.14159265358979328
 print*,b
 b=3.14159265358979328d0
 print*,b
 c=3.14159265358979328q0
 print*,c

 end




 
 
