module str2int_mod ! from https://stackoverflow.com/questions/24071722/
    contains 
          elemental subroutine str2int(str,int,stat)
          implicit none
          ! Arguments
          character(len=*),intent(in) :: str
          integer,intent(out)         :: int
          integer,intent(out)         :: stat
    
          read(str,*,iostat=stat)  int
          end subroutine str2int
end module

PROGRAM str2int_test
    use str2int_mod
    !     IMPLICIT NONE
    !     .. Parameters ..
    !     getting allocated memory from command line
    INTEGER i, status
    CHARACTER*32 arg

    CALL getarg(1, arg)
    CALL str2int(arg, i, status)

    if ( status == 0 ) then
        print *,i
    else
        print *,'Conversion of string ',i,' failed!'
    endif
END PROGRAM str2int_test