program ifendif

    implicit none

    ! Declare the variables, ask for input
    real :: x, y
    integer :: choice

    print *, "What are your x and y values?"
    read *, x
    read *, y 

    ! Ask user to choose from a list of options
    print *, "Please choose an option: "
    print *, "1. Addition"
    print *, "2. Subtraction"
    print *, "3. Multiplication"
    print *, "4. Division"

    read *, choice 

    ! Act accordingly
    if (choice == 1) print *, "Your result is", x + y
    if (choice == 2) print *, "Your result is", x - y
    if (choice == 3) print *, "Your result is", x * y
    if (choice == 4) print *, "Your result is", x / y ! The guide suggested putting a tertiary variable here, which is silly imo

    ! If I wanted to do this with more complicated if statements, the syntax is:
    ! if (choice == 5) then
    !     print *, "This is a more complicated if statement involving more than 1 line."
    !     print *, "Bow before its might and wonder."
    ! end if

    ! Additionally: if I want to compare a value to 0, like (x == 0), FORTRAN can fail to evaluate it to a hard 0. Thus, I should do something like:
    ! if (abs(x) < 0.000001) then ... ! where abs() takes the absolute value of something.

end program ifendif