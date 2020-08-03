program writedata 

    implicit none

    ! Declare vars
    integer :: x, y, z

    print *, "What are your three vars for output into a file?"
    read *, x
    read *, y 
    read *, z

    ! Main
    open(10, file="output.txt")
    write(10, *) x, y, z

end program writedata