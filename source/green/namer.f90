program namer

    implicit none
    
    ! Declare the variables
    Character :: name * 10

    ! Ask me to write my name.
    print *, "What is your name? "
    
    ! Recognize the typing from the terminal.
    read *, name
    
    ! Then, print my name.
    print *, "My name is ", name

end program namer