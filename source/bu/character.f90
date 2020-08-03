
 implicit none

 integer       :: n 
 character     :: a,b
 character(10) :: c

 a='A'
 b='B'
 c=a//b//b//a
 n=len_trim(c)

 print*,'Number of characters in c',n
 print*,c(1:n),' ',c(2:3),' ',c(1:1)
 print*,iachar(a),iachar(b)
 print*,char(7),char(7),char(67),char(68)

 end

 
