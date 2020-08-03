
 implicit none

 integer :: int

 print*,'Give an integer between 1 and 99'; read*,int
 select case (int)
 case (:0, 100:)
   print*,'Read the instructions more carefully next time! Good bye.'
 case (8,88)
   print*,'A lucky number; Congratulations!'
 case (13)
   print*,'Bad luck...not a good number; beware!'
 case default
   print*,'Nothing special with this number,'
   select case (mod(int,2)==0)
     case (.true.)
       print*,'but I can tell you that it is an even number.'
     case (.false.)
       print*,'but I can tell you that it is an odd number.'
   end select
 end select

 end



