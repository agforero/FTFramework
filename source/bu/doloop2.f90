
 implicit none

 integer :: i,n

 print*,'Give highest number n to be squared'; read*,n

 i=0
 do
   i=i+1
   print*,i**2
   IF (i.EQ.n) exit 
 end do

 end

