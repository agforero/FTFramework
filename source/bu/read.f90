 implicit none

 integer :: nlines
 character     :: c
 character(16) :: fname

 print*,'Give a file:';read*,fname
 open(1,file=fname,status='old')
 nlines=0
 do
   read(1,*,end=10)c; print*,c
   nlines=nlines+1
 end do
 10 close(1)
 print*; print*,'The number of non-empty lines in the file is:',nlines

 end
