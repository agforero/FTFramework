program helloworld
integer :: i
character(len=32) :: arg

call getarg(1, arg)
print *, arg

do i = 1, iargc()
        call getarg(i, arg)
        print *, arg
end do

end program helloworld
