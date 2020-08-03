program writedatamulti

    implicit none

    integer :: x, y, z, a, b, c

    print *, "What are your x y and z vals?"
    read *, x
    read *, y 
    read *, z

    ! Main 
    open(10, file="outputmulti.txt")
    
    do a=0, x
        do b=0, y 
            do c=0, z
                write(10,*), a, b, c
            end do
        end do
    end do

end program writedatamulti