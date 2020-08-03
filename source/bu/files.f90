
 implicit none

 integer i,nfiles
 character(10) fname

 10 print*,'Number of files to open'; read*,nfiles
 if (nfiles > 99) then
   print*,"That's a little excessive; try something reasonable!"
   goto 10
 end if
 fname='file00.txt'
 do i=1,nfiles
   fname(6:6)=achar(48+mod(i,10))
   fname(5:5)=achar(48+i/10)
   open(1,file=fname,status='new')
   write(1,*)'This is file # ',i
   close(1)
 enddo

 end

