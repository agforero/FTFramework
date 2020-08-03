
 implicit none

 character(60) :: fname,line

 write(*,'(a)',advance='no')'Give a file: ';read*,fname
 print*
 open(1,file=fname,status='old')
 do 
   read(1,*,end=10)line; print*,line
 end do
 10 close(1)
 open(1,file=fname,status='old')
 do 
   read(1,'(a)',end=20)line; print*,line
 end do
 20 close(1)

 end




