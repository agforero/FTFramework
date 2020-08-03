
 implicit none

 integer :: i,a,b,c,p
 character(32) :: bits

 print*,'Give integer a'; read*,a
 print*,'Give integer b'; read*,b
 print*,'Bit position p for single-bit functions on a'; read*,p

 print*,'a          : ',bits(a),'  =  ',a
 print*,'b          : ',bits(b),'  =  ',b
 c=ibset(a,p); print*,'ibset(a,p) : ',bits(c),'  =  ',c
 c=ibclr(a,p); print*,'ibclr(a,p) : ',bits(c),'  =  ',c
 c=ishft(a,p); print*,'ishft(a,p) : ',bits(c),'  =  ',c
 c=iand(a,b);  print*,'iand(a,b)  : ',bits(c),'  =  ',c
 c=ior(a,b);   print*,'ior(a,b)   : ',bits(c),'  =  ',c
 c=ieor(a,b);  print*,'ieor(a,b)  : ',bits(c),'  =  ',c

 end

 function bits(int)

 implicit none

 integer :: i,int
 character(32) :: bits

 bits='00000000000000000000000000000000'
 do i=0,31
   if (btest(int,i)) bits(32-i:32-i)='1'
 end do

 end function bits
