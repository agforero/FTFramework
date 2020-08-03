program divide

    implicit none

    ! Create vars
    real :: x, y, ans

    ! Ask for em
    print *, "What are your x and y vals?"

    ! Assign
    read *, x
    read *, y 

    ans = x / y 

    ! Display
    print *, "Your result is ", ans

end program divide