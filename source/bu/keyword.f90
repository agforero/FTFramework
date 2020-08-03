 module test

 contains

    subroutine keywordsub(a,b)
    integer, optional :: a
    integer, optional :: b
    if (present(a)) write(*,*)'a = ',a
    if (present(b)) write(*,*)'b = ',b
    end subroutine keywordsub

 end module test

program testkeyword
use test

integer :: arg1,arg2

write(*,*)'Give 2 integers:'
read(*,*)arg1,arg2

write(*,*)
call keywordsub(arg1)
write(*,*)
call keywordsub(b=arg2)
write(*,*)
call keywordsub(a=arg1,b=arg2)

end program testkeyword
