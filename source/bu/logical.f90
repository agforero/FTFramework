 
 implicit none
 
 logical :: a,b

 print*,'Give values (T/F) for a and b'
 read*,a,b
 print*,a.or.b,a.and.b,a.eqv.b,a.neqv.b,.not.a

 end
