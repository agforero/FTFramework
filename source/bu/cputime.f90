
 implicit none

 integer :: i,nloop
 real(8) :: sum
 real    :: time0,time1
 
 print*,'Number of operations in each loop'; read*,nloop

 sum=0.0_8; call cpu_time(time0)
 do i=1,nloop
   sum=sum+dble(i)*dble(i)
 end do
 print*,sum
 call cpu_time(time1); print*,'Time used for s=s+i*i     : ',time1-time0
 sum=0.0_8; call cpu_time(time0)
 do i=1,nloop
   sum=sum+sqrt(dble(i))
 end do
 print*,sum
 call cpu_time(time1); print*,'Time used for s=s+sqrt(i) : ',time1-time0

 sum=0.0_8; call cpu_time(time0)
 do i=1,nloop
   sum=sum+log(dble(i))
 end do
 print*,sum
 call cpu_time(time1); print*,'Time used for s=s+log(i)  : ',time1-time0

 sum=0.0_8; call cpu_time(time0)
 do i=1,nloop
   sum=sum+exp(1.d0/dble(i))
 end do
 print*,sum
 call cpu_time(time1); print*,'Time used for s=s+exp(1/i)  : ',time1-time0

 sum=0.0_8; call cpu_time(time0)
 do i=1,nloop
   sum=sum+cos(dble(i))
 end do
 print*,sum
 call cpu_time(time1); print*,'Time used for s=s+cos(i)  : ',time1-time0

 end


 
