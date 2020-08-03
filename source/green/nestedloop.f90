program nestedloop

    implicit none

    integer :: x, y

    do x = 0, 5
        do y = 0, 5
            print *, x, y
        end do
    end do ! hmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

end program nestedloop