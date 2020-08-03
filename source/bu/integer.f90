
 implicit none

 integer :: i

 i=(-2)**31  ! the 31 power if -2 is OK
! i=-2**31   ! - (the 31 power of 2) is not allowed by the compiler
 print*,i
 i=i-1
 print*,i

 end
