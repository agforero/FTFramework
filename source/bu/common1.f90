
 implicit none

 integer :: a,b
 common/block_1/a,b

 integer c,n

 print*,'Give integers a and b'; read*,a,b
 print*,'How many strange operations should be performed?'; read*,n

 call strangeoperation(n,c)
 print*,'The strange operations have resulted in the value ',c

 end

 subroutine strangeoperation(n,strange)

 implicit none

 integer :: a,b
 common/block_1/a,b

 integer i,n,strange

 strange=0
 do i=1,n
   strange=strange+(strange-a)**2+a*b*strange
   if (abs(strange) > 100000) strange=strange/1000+10*i
 enddo

 end subroutine strangeoperation



