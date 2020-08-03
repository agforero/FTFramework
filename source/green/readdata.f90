program readdata

    implicit none

    ! Declaring vars
    integer :: x, y, z

    ! Main
    open(10, file="output.txt")
    read(10, *) x, y, z

    print *, "The numbers contained in the file are", x, y, z

end program readdata