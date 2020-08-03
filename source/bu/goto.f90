
 implicit none

 integer :: i,i2,sqmax

 print*,'Give upper limit for squares to be produced'; read*,sqmax

 i=0
 10 i=i+1
   i2=i**2
   if (i2 < sqmax) then
     print*,'Integer: ',i,'  Its square: ',i2
     go to 10
   end if
 
 end

    
 