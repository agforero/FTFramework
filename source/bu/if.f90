
 implicit none

 integer :: int

 print*,'Give an integer between 1 and 99'; read*,int
 if (int<1.or.int>99) then
   print*,'Read the instructions more carefully next time! Good bye.'
 else if (int==8.or.int==88) then
   print*,'A lucky number; Congratulations!'
 else if (int==13.or.int==4) then
   print*,'Bad luck...not a good number; beware!'
 else
   print*,'Nothing special with this number, '
   if (mod(int,2)==0) then
     print*,'but I can tell you that it is an even number.'
   else
     print*,'but I can tell you that it is an odd number.'
   end if
 end if

 end
