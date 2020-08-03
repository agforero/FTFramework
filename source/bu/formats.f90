
 implicit none

 integer :: a(4)
 real(8) :: b(4)

 a(1)=1; a(2)=10; a(3)=100; a(4)=1000
 b(1)=1.d0; b(2)=1.d1; b(3)=1.d2; b(4)=1.d3

 print'(4i5)',a
 write(*,'(4i5)'),a
 write(*,10)a
 10 format(4i5)
 print'(4i3)',a
 print'(a,i1,a,i2,a,i3)',' one:',a(1),' ten:',a(2),' hundred:',a(3)
 print'(6f12.6)',b


 end
